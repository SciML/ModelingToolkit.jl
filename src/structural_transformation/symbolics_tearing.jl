# N.B. assumes `slist` and `dlist` are unique
function substitution_graph(graph, slist, dlist, var_eq_matching)
    ns = length(slist)
    nd = length(dlist)
    ns == nd || error("internal error")
    newgraph = BipartiteGraph(ns, nd)
    erename = uneven_invmap(nsrcs(graph), slist)
    vrename = uneven_invmap(ndsts(graph), dlist)
    for e in 𝑠vertices(graph)
        ie = erename[e]
        ie == 0 && continue
        for v in 𝑠neighbors(graph, e)
            iv = vrename[v]
            iv == 0 && continue
            add_edge!(newgraph, ie, iv)
        end
    end

    newmatching = Matching(ns)
    for (v, e) in enumerate(var_eq_matching)
        (e === unassigned || e === SelectedState()) && continue
        iv = vrename[v]
        ie = erename[e]
        iv == 0 && continue
        ie == 0 && error("internal error")
        newmatching[iv] = ie
    end

    return DiCMOBiGraph{true}(newgraph, complete(newmatching))
end

function var_derivative_graph!(s::SystemStructure, v::Int)
    sg = g = add_vertex!(s.graph, DST)
    var_diff = add_vertex!(s.var_to_diff)
    add_edge!(s.var_to_diff, v, var_diff)
    s.solvable_graph === nothing || (sg = add_vertex!(s.solvable_graph, DST))
    @assert sg == g == var_diff
    return var_diff
end

function var_derivative!(ts::TearingState{ODESystem}, v::Int)
    s = ts.structure
    var_diff = var_derivative_graph!(s, v)
    sys = ts.sys
    D = Differential(get_iv(sys))
    push!(ts.fullvars, D(ts.fullvars[v]))
    return var_diff
end

function eq_derivative_graph!(s::SystemStructure, eq::Int)
    add_vertex!(s.graph, SRC)
    s.solvable_graph === nothing || add_vertex!(s.solvable_graph, SRC)
    # the new equation is created by differentiating `eq`
    eq_diff = add_vertex!(s.eq_to_diff)
    add_edge!(s.eq_to_diff, eq, eq_diff)
    return eq_diff
end

function eq_derivative!(ts::TearingState{ODESystem}, ieq::Int)
    s = ts.structure

    eq_diff = eq_derivative_graph!(s, ieq)

    sys = ts.sys
    D = Differential(get_iv(sys))
    eq = equations(ts)[ieq]
    eq = ModelingToolkit.expand_derivatives(0 ~ D(eq.rhs - eq.lhs))
    push!(equations(ts), eq)
    # Analyze the new equation and update the graph/solvable_graph
    # First, copy the previous incidence and add the derivative terms.
    # That's a superset of all possible occurrences. find_solvables! will
    # remove those that doesn't actually occur.
    eq_diff = length(equations(ts))
    for var in 𝑠neighbors(s.graph, ieq)
        add_edge!(s.graph, eq_diff, var)
        add_edge!(s.graph, eq_diff, s.var_to_diff[var])
    end
    s.solvable_graph === nothing ||
        find_eq_solvables!(ts, eq_diff; may_be_zero = true, allow_symbolic = true)

    return eq_diff
end

function tearing_sub(expr, dict, s)
    expr = ModelingToolkit.fixpoint_sub(expr, dict)
    s ? simplify(expr) : expr
end

function full_equations(sys::AbstractSystem; simplify = false)
    empty_substitutions(sys) && return equations(sys)
    substitutions = get_substitutions(sys)
    substitutions.subed_eqs === nothing || return substitutions.subed_eqs
    @unpack subs = substitutions
    solved = Dict(eq.lhs => eq.rhs for eq in subs)
    neweqs = map(equations(sys)) do eq
        if isdiffeq(eq)
            return tearing_sub(eq.lhs, solved, simplify) ~ tearing_sub(eq.rhs, solved,
                                                                       simplify)
        else
            if !(eq.lhs isa Number && eq.lhs == 0)
                eq = 0 ~ eq.rhs - eq.lhs
            end
            rhs = tearing_sub(eq.rhs, solved, simplify)
            if rhs isa Symbolic
                return 0 ~ rhs
            else # a number
                error("tearing failled because the system is singular")
            end
        end
        eq
    end
    substitutions.subed_eqs = neweqs
    return neweqs
end

function tearing_substitution(sys::AbstractSystem; kwargs...)
    neweqs = full_equations(sys::AbstractSystem; kwargs...)
    @set! sys.eqs = neweqs
    @set! sys.substitutions = nothing
end

function tearing_assignments(sys::AbstractSystem)
    if empty_substitutions(sys)
        assignments = []
        deps = Int[]
        sol_states = Code.LazyState()
    else
        @unpack subs, deps = get_substitutions(sys)
        assignments = [Assignment(eq.lhs, eq.rhs) for eq in subs]
        sol_states = Code.NameState(Dict(eq.lhs => Symbol(eq.lhs) for eq in subs))
    end
    return assignments, deps, sol_states
end

function solve_equation(eq, var, simplify)
    rhs = value(solve_for(eq, var; simplify = simplify, check = false))
    occursin(var, rhs) && throw(EquationSolveErrors(eq, var, rhs))
    var ~ rhs
end

function substitute_vars!(graph::BipartiteGraph, subs, cache = Int[], callback! = nothing;
                          exclude = ())
    for su in subs
        su === nothing && continue
        v, v′ = su
        eqs = 𝑑neighbors(graph, v)
        # Note that the iterator is not robust under deletion and
        # insertion. Hence, we have a copy here.
        resize!(cache, length(eqs))
        for eq in copyto!(cache, eqs)
            eq in exclude && continue
            rem_edge!(graph, eq, v)
            add_edge!(graph, eq, v′)
            callback! !== nothing && callback!(eq, su)
        end
    end
    graph
end

function to_mass_matrix_form(neweqs, ieq, graph, fullvars, isdervar::F,
                             var_to_diff) where {F}
    eq = neweqs[ieq]
    if !(eq.lhs isa Number && eq.lhs == 0)
        eq = 0 ~ eq.rhs - eq.lhs
    end
    rhs = eq.rhs
    if rhs isa Symbolic
        # Check if the RHS is solvable in all state derivatives and if those
        # the linear terms for them are all zero. If so, move them to the
        # LHS.
        dervar::Union{Nothing, Int} = nothing
        for var in 𝑠neighbors(graph, ieq)
            if isdervar(var)
                if dervar !== nothing
                    error("$eq has more than one differentiated variable!")
                end
                dervar = var
            end
        end
        dervar === nothing && return (0 ~ rhs), dervar
        new_lhs = var = fullvars[dervar]
        # 0 ~ a * D(x) + b
        # D(x) ~ -b/a
        a, b, islinear = linear_expansion(rhs, var)
        if !islinear
            return (0 ~ rhs), nothing
        end
        new_rhs = -b / a
        return (new_lhs ~ new_rhs), invview(var_to_diff)[dervar]
    else # a number
        if abs(rhs) > 100eps(float(rhs))
            @warn "The equation $eq is not consistent. It simplifed to 0 == $rhs."
        end
        return nothing
    end
end

function tearing_reassemble(state::TearingState, var_eq_matching; simplify = false)
    @unpack fullvars, sys = state
    @unpack solvable_graph, var_to_diff, eq_to_diff, graph = state.structure

    neweqs = collect(equations(state))
    # substitution utilities
    idx_buffer = Int[]
    sub_callback! = let eqs = neweqs, fullvars = fullvars
        (ieq, s) -> begin
            neweq = substitute(eqs[ieq], fullvars[s[1]] => fullvars[s[2]])
            eqs[ieq] = neweq
        end
    end

    # Terminology and Definition:
    #
    # A general DAE is in the form of `F(u'(t), u(t), p, t) == 0`. We can
    # characterize variables in `u(t)` into two classes: differential variables
    # (denoted `v(t)`) and algebraic variables (denoted `z(t)`). Differential
    # variables are marked as `SelectedState` and they are differentiated in the
    # DAE system, i.e. `v'(t)` are all the variables in `u'(t)` that actually
    # appear in the system. Algebraic variables are variables that are not
    # differential variables.
    #
    # Dummy derivatives may determine that some differential variables are
    # algebraic variables in disguise. The derivative of such variables are
    # called dummy derivatives.

    # Step 1:
    # Replace derivatives of non-selected states by dummy derivatives

    possible_x_t = Dict()
    oldobs = observed(sys)
    for (i, eq) in enumerate(oldobs)
        lhs = eq.lhs
        rhs = eq.rhs
        isdifferential(lhs) && continue
        # TODO: should we hanlde negative alias as well?
        isdifferential(rhs) || continue
        possible_x_t[rhs] = i, lhs
    end

    #removed_eqs = Int[]
    #removed_vars = Int[]
    removed_obs = Int[]
    diff_to_var = invview(var_to_diff)
    for var in 1:length(fullvars)
        dv = var_to_diff[var]
        dv === nothing && continue
        if var_eq_matching[var] !== SelectedState()
            dd = fullvars[dv]
            if (i_v_t = get(possible_x_t, dd, nothing)) === nothing
                v_t = diff2term(unwrap(dd))
            else
                idx, v_t = i_v_t
                push!(removed_obs, idx)
            end
            for eq in 𝑑neighbors(graph, dv)
                neweqs[eq] = substitute(neweqs[eq], fullvars[dv] => v_t)
            end
            fullvars[dv] = v_t
            # update the structural information
            diff_to_var[dv] = nothing
        end
    end

    # `SelectedState` information is no longer needed past here. State selection
    # is done. All non-differentiated variables are algebraic variables, and all
    # variables that appear differentiated are differential variables.

    ### extract partition information
    is_solvable = let solvable_graph = solvable_graph
        (eq, iv) -> eq isa Int && iv isa Int && BipartiteEdge(eq, iv) in solvable_graph
    end

    # if var is like D(x)
    isdervar = let diff_to_var = diff_to_var
        var -> diff_to_var[var] !== nothing
    end

    # There are three cases where we want to generate new variables to convert
    # the system into first order (semi-implicit) ODEs.
    #
    # 1. To first order:
    # Whenever higher order differentiated variable like `D(D(D(x)))` appears,
    # we introduce new variables `x_t`, `x_tt`, and `x_ttt` and new equations
    # ```
    # D(x_tt) = x_ttt
    # D(x_t) = x_tt
    # D(x) = x_t
    # ```
    # and replace `D(x)` to `x_t`, `D(D(x))` to `x_tt`, and `D(D(D(x)))` to
    # `x_ttt`.
    #
    # 2. To implicit to semi-implicit ODEs:
    # 2.1: Unsolvable derivative:
    # If one derivative variable `D(x)` is unsolvable in all the equations it
    # appears in, then we introduce a new variable `x_t`, a new equation
    # ```
    # D(x) ~ x_t
    # ```
    # and replace all other `D(x)` to `x_t`.
    #
    # 2.2: Solvable derivative:
    # If one derivative variable `D(x)` is solvable in at least one of the
    # equations it appears in, then we introduce a new variable `x_t`. One of
    # the solvable equations must be in the form of `0 ~ L(D(x), u...)` and
    # there exists a function `l` such that `D(x) ~ l(u...)`. We should replace
    # it to
    # ```
    # 0 ~ x_t - l(u...)
    # D(x) ~ x_t
    # ```
    # and replace all other `D(x)` to `x_t`.
    #
    # Observe that we don't need to actually introduce a new variable `x_t`, as
    # the above equations can be lowered to
    # ```
    # x_t := l(u...)
    # D(x) ~ x_t
    # ```
    # where `:=` denotes assignment.
    #
    # As a final note, in all the above cases where we need to introduce new
    # variables and equations, don't add them when they already exist.

    if ModelingToolkit.has_iv(state.sys)
        iv = get_iv(state.sys)
        D = Differential(iv)
    else
        iv = D = nothing
    end
    nvars = ndsts(graph)
    processed = falses(nvars)
    subinfo = NTuple{3, Int}[]
    for i in 1:nvars
        processed[i] && continue

        v = i
        # descend to the bottom of differentiation chain
        while diff_to_var[v] !== nothing
            v = diff_to_var[v]
        end

        # `v` is now not differentiated at level 0.
        diffvar = v
        processed[v] = true
        level = 0
        order = 0
        isimplicit = false
        # ascend to the top of differentiation chain
        while true
            eqs_with_v = 𝑑neighbors(graph, v)
            if !isempty(eqs_with_v)
                order = level
                isimplicit = length(eqs_with_v) > 1 || !is_solvable(only(eqs_with_v), v)
            end
            if v <= length(processed)
                processed[v] = true
            end
            var_to_diff[v] === nothing && break
            v = var_to_diff[v]
            level += 1
        end

        # `diffvar` is a order `order` variable
        (isimplicit || order > 1) || continue

        # add `D(t) ~ x_t` etc
        subs = Dict()
        ogx = x = fullvars[diffvar] # x
        ogidx = xidx = diffvar
        # We shouldn't apply substitution to `order_lowering_eqs`
        order_lowering_eqs = BitSet()
        for o in 1:order
            # D(x) ~ x_t
            ogidx = var_to_diff[ogidx]

            dx_idx = var_to_diff[xidx]
            if dx_idx === nothing
                dx = D(x)
                push!(fullvars, dx)
                dx_idx = add_vertex!(var_to_diff)
                add_vertex!(graph, DST)
                add_vertex!(solvable_graph, DST)
                @assert dx_idx == ndsts(graph) == length(fullvars)
                push!(var_eq_matching, unassigned)

                var_to_diff[xidx] = dx_idx
            else
                dx = fullvars[dx_idx]
            end

            if (i_x_t = get(possible_x_t, dx, nothing)) === nothing &&
               (ogidx !== nothing &&
                (i_x_t = get(possible_x_t, fullvars[ogidx], nothing)) === nothing)
                x_t = ModelingToolkit.lower_varname(ogx, iv, o)
            else
                idx, x_t = i_x_t
                push!(removed_obs, idx)
            end
            push!(fullvars, x_t)
            x_t_idx = add_vertex!(var_to_diff)
            add_vertex!(graph, DST)
            add_vertex!(solvable_graph, DST)
            @assert x_t_idx == ndsts(graph) == length(fullvars)
            push!(var_eq_matching, unassigned)

            push!(neweqs, dx ~ x_t)
            eq_idx = add_vertex!(eq_to_diff)
            push!(order_lowering_eqs, eq_idx)
            add_vertex!(graph, SRC)
            add_vertex!(solvable_graph, SRC)
            @assert eq_idx == nsrcs(graph) == length(neweqs)

            add_edge!(solvable_graph, eq_idx, x_t_idx)
            add_edge!(solvable_graph, eq_idx, dx_idx)
            add_edge!(graph, eq_idx, x_t_idx)
            add_edge!(graph, eq_idx, dx_idx)
            # We use this info to substitute all `D(D(x))` or `D(x_t)` except
            # the `D(D(x)) ~ x_tt` equation to `x_tt`.
            #              D(D(x))  D(x_t)    x_tt
            push!(subinfo, (ogidx, dx_idx, x_t_idx))

            # D(x_t) ~ x_tt
            x = x_t
            xidx = x_t_idx
        end

        # Go backward from high order to lower order so that we substitute
        # something like `D(D(x)) -> x_tt` first, otherwise we get `D(x_t)`
        # which would be hard to fix up before we finish lower the order of
        # variable `x`.
        for (ogidx, dx_idx, x_t_idx) in Iterators.reverse(subinfo)
            # We need a loop here because both `D(D(x))` and `D(x_t)` need to be
            # substituted to `x_tt`.
            for idx in (ogidx, dx_idx)
                subidx = ((idx => x_t_idx),)
                # This handles case 2.2
                if var_eq_matching[idx] isa Int
                    var_eq_matching[x_t_idx] = var_eq_matching[idx]
                end
                substitute_vars!(graph, subidx, idx_buffer, sub_callback!;
                                 exclude = order_lowering_eqs)
                substitute_vars!(solvable_graph, subidx, idx_buffer;
                                 exclude = order_lowering_eqs)
            end
        end
        empty!(subinfo)
        empty!(subs)
    end

    # Will reorder equations and states to be:
    # [diffeqs; ...]
    # [diffvars; ...]
    # such that the mass matrix is:
    # [I  0
    #  0  0].
    diffeq_idxs = Int[]
    algeeq_idxs = Int[]
    diff_eqs = Equation[]
    alge_eqs = Equation[]
    diff_vars = Int[]
    subeqs = Equation[]
    solved_equations = Int[]
    solved_variables = Int[]
    # Solve solvable equations
    neqs = nsrcs(graph)
    for (ieq, iv) in enumerate(invview(var_eq_matching))
        ieq > neqs && break
        if is_solvable(ieq, iv)
            # We don't solve differential equations, but we will need to try to
            # convert it into the mass matrix form.
            # We cannot solve the differential variable like D(x)
            if isdervar(iv)
                # TODO: what if `to_mass_matrix_form(ieq)` returns `nothing`?
                eq, diffidx = to_mass_matrix_form(neweqs, ieq, graph, fullvars, isdervar,
                                                  var_to_diff)
                push!(diff_eqs, eq)
                push!(diffeq_idxs, ieq)
                push!(diff_vars, diffidx)
                continue
            end
            eq = neweqs[ieq]
            var = fullvars[iv]
            residual = eq.lhs - eq.rhs
            a, b, islinear = linear_expansion(residual, var)
            @assert islinear
            # 0 ~ a * var + b
            # var ~ -b/a
            if ModelingToolkit._iszero(a)
                @warn "Tearing: $eq is a singular equation!"
                #push!(removed_eqs, ieq)
                #push!(removed_vars, iv)
            else
                rhs = -b / a
                neweq = var ~ simplify ? Symbolics.simplify(rhs) : rhs
                push!(subeqs, neweq)
                push!(solved_equations, ieq)
                push!(solved_variables, iv)
            end
        else
            eq, diffidx = to_mass_matrix_form(neweqs, ieq, graph, fullvars, isdervar,
                                              var_to_diff)
            if diffidx === nothing
                push!(alge_eqs, eq)
                push!(algeeq_idxs, ieq)
            else
                push!(diff_eqs, eq)
                push!(diffeq_idxs, ieq)
                push!(diff_vars, diffidx)
            end
        end
    end
    # TODO: BLT sorting
    neweqs = [diff_eqs; alge_eqs]
    inveqsperm = [diffeq_idxs; algeeq_idxs]
    eqsperm = zeros(Int, nsrcs(graph))
    for (i, v) in enumerate(inveqsperm)
        eqsperm[v] = i
    end
    diff_vars_set = BitSet(diff_vars)
    if length(diff_vars_set) != length(diff_vars)
        error("Tearing internal error: lowering DAE into semi-implicit ODE failed!")
    end
    invvarsperm = [diff_vars;
                   setdiff!(setdiff(1:ndsts(graph), diff_vars_set),
                            BitSet(solved_variables))]
    varsperm = zeros(Int, ndsts(graph))
    for (i, v) in enumerate(invvarsperm)
        varsperm[v] = i
    end

    if isempty(solved_equations)
        deps = Vector{Int}[]
    else
        subgraph = substitution_graph(graph, solved_equations, solved_variables,
                                      var_eq_matching)
        toporder = topological_sort_by_dfs(subgraph)
        subeqs = subeqs[toporder]
        # Find the dependency of solved variables. We will need this for ODAEProblem
        invtoporder = invperm(toporder)
        deps = [Int[invtoporder[n]
                    for n in neighborhood(subgraph, j, Inf, dir = :in) if n != j]
                for j in toporder]
    end

    # Contract the vertices in the structure graph to make the structure match
    # the new reality of the system we've just created.
    graph = contract_variables(graph, var_eq_matching, varsperm, eqsperm,
                               length(solved_variables))

    # Update system
    new_var_to_diff = complete(DiffGraph(length(invvarsperm)))
    for (v, d) in enumerate(var_to_diff)
        v′ = varsperm[v]
        (v′ > 0 && d !== nothing) || continue
        d′ = varsperm[d]
        new_var_to_diff[v′] = d′ > 0 ? d′ : nothing
    end
    new_eq_to_diff = complete(DiffGraph(length(inveqsperm)))
    for (v, d) in enumerate(eq_to_diff)
        v′ = eqsperm[v]
        (v′ > 0 && d !== nothing) || continue
        d′ = eqsperm[d]
        new_eq_to_diff[v′] = d′ > 0 ? d′ : nothing
    end

    var_to_diff = new_var_to_diff
    eq_to_diff = new_eq_to_diff
    diff_to_var = invview(var_to_diff)

    @set! state.structure.graph = complete(graph)
    @set! state.structure.var_to_diff = var_to_diff
    @set! state.structure.eq_to_diff = eq_to_diff
    @set! state.fullvars = fullvars = fullvars[invvarsperm]

    sys = state.sys
    @set! sys.eqs = neweqs
    @set! sys.states = [v for (i, v) in enumerate(fullvars) if diff_to_var[i] === nothing]
    removed_obs_set = BitSet(removed_obs)
    var_to_idx = Dict(reverse(en) for en in enumerate(fullvars))
    # Make sure differentiated variables don't appear in observed equations
    for (dx, (idx, lhs)) in possible_x_t
        idx in removed_obs_set && continue
        # Because it's a differential variable, and by sorting, its
        # corresponding differential equation would have the same index.
        dxidx = get(var_to_idx, dx, nothing)
        # TODO: use alias graph to handle the dxidx === nothing case for
        # mechanical systems.
        dxidx === nothing && continue
        eqidx = diff_to_var[dxidx]
        oldobs[idx] = (lhs ~ neweqs[eqidx].rhs)
    end
    deleteat!(oldobs, sort!(removed_obs))
    @set! sys.observed = [oldobs; subeqs]
    @set! sys.substitutions = Substitutions(subeqs, deps)
    @set! state.sys = sys
    @set! sys.tearing_state = state

    return invalidate_cache!(sys)
end

function tearing(state::TearingState; kwargs...)
    state.structure.solvable_graph === nothing && find_solvables!(state; kwargs...)
    complete!(state.structure)
    @unpack graph = state.structure
    algvars = BitSet(findall(v -> isalgvar(state.structure, v), 1:ndsts(graph)))
    aeqs = algeqs(state.structure)
    var_eq_matching′ = tear_graph_modia(state.structure;
                                        varfilter = var -> var in algvars,
                                        eqfilter = eq -> eq in aeqs)
    var_eq_matching = Matching{Union{Unassigned, SelectedState}}(var_eq_matching′)
    for var in 1:ndsts(graph)
        if isdiffvar(state.structure, var)
            var_eq_matching[var] = SelectedState()
        end
    end
    var_eq_matching
end

"""
    tearing(sys; simplify=false)

Tear the nonlinear equations in system. When `simplify=true`, we simplify the
new residual residual equations after tearing. End users are encouraged to call [`structural_simplify`](@ref)
instead, which calls this function internally.
"""
function tearing(sys::AbstractSystem; simplify = false)
    state = TearingState(sys)
    var_eq_matching = tearing(state)
    invalidate_cache!(tearing_reassemble(state, var_eq_matching; simplify = simplify))
end

"""
    partial_state_selection(sys; simplify=false)

Perform partial state selection and tearing.
"""
function partial_state_selection(sys; simplify = false)
    state = TearingState(sys)
    var_eq_matching = partial_state_selection_graph!(state)

    tearing_reassemble(state, var_eq_matching; simplify = simplify)
end

"""
    dummy_derivative(sys)

Perform index reduction and use the dummy derivative technique to ensure that
the system is balanced.
"""
function dummy_derivative(sys, state = TearingState(sys); simplify = false, kwargs...)
    function jac(eqs, vars)
        symeqs = EquationsView(state)[eqs]
        Symbolics.jacobian((x -> x.rhs).(symeqs), state.fullvars[vars])
    end
    var_eq_matching = dummy_derivative_graph!(state, jac; kwargs...)
    tearing_reassemble(state, var_eq_matching; simplify = simplify)
end
