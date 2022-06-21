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

function var_derivative!(ts::TearingState{ODESystem}, v::Int)
    sys = ts.sys
    s = ts.structure
    D = Differential(get_iv(sys))
    s.solvable_graph === nothing || add_vertex!(s.solvable_graph, DST)
    push!(ts.fullvars, D(ts.fullvars[v]))
end

function eq_derivative!(ts::TearingState{ODESystem}, ieq::Int)
    sys = ts.sys
    s = ts.structure
    D = Differential(get_iv(sys))
    eq = equations(ts)[ieq]
    eq = ModelingToolkit.expand_derivatives(0 ~ D(eq.rhs - eq.lhs))
    s.solvable_graph === nothing || add_vertex!(s.solvable_graph, SRC)
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

# From the index of `D(x)` find the equation `D(x) ~ x_t` and the variable
# `x_t`.
function has_order_lowering_eq_var(eqs, fullvars, graph, var_to_diff, dx_idx)::Union{Nothing, NTuple{2, Int}}
    diff_to_var = invview(var_to_diff)
    diff_to_var[dx_idx] === nothing && return nothing

    dx = fullvars[dx_idx]
    for eq in 𝑑neighbors(graph, dx_idx)
        vs = 𝑠neighbors(graph, eq)
        length(vs) == 2 || continue
        maybe_x_t_idx = vs[1] == dx_idx ? vs[2] : vs[1]
        # TODO: should we follow the differentiation chain? I.e. recurse until
        # all reachable variables are explored or `diff_to_var[maybe_x_t_idx] === nothing`
        diff_to_var[maybe_x_t_idx] === nothing || continue
        maybe_x_t = fullvars[maybe_x_t_idx]
        difference = (eqs[eq].lhs - eqs[eq].rhs) - (dx - maybe_x_t)
        # if `eq` is in the form of `D(x) ~ x_t`
        if ModelingToolkit._iszero(difference)
            # TODO: reduce systems with multiple order lowering `eq` and `var`
            # as well.
            return eq, maybe_x_t_idx
        end
    end
    return nothing
end

function var2var_t_map(state::TearingState)
    fullvars = state.fullvars
    @unpack var_to_diff, graph = state.structure
    eqs = equations(state)
    @info "" eqs
    var2var_t = Vector{Union{Nothing, NTuple{2, Int}}}(undef, ndsts(graph))
    for v in 1:ndsts(graph)
        var2var_t[v] = has_order_lowering_eq_var(eqs, fullvars, graph, var_to_diff, v)
    end
    var2var_t
end

function substitute_vars!(graph::BipartiteGraph, subs, cache=Int[], callback! = nothing; exclude = ())
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

function tearing_reassemble(state::TearingState, var_eq_matching, var2var_t = var2var_t_map(state); simplify = false)
    fullvars = state.fullvars
    @unpack solvable_graph, var_to_diff, eq_to_diff, graph = state.structure

    neweqs = collect(equations(state))
    # substitution utilities
    idx_buffer = Int[]
    sub_callback! = let eqs = neweqs, fullvars = fullvars
        (ieq, s) -> eqs[ieq] = substitute(eqs[ieq], fullvars[s[1]] => fullvars[s[2]])
    end

    @info "" neweqs

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

    remove_eqs = Int[]
    diff_to_var = invview(var_to_diff)
    for var in 1:length(fullvars)
        dv = var_to_diff[var]
        dv === nothing && continue
        if var_eq_matching[var] !== SelectedState()
            dd = fullvars[dv]
            # convert `D(x)` to `x_t` (don't rely on the specific spelling of
            # the name)
            eq_var_t = var2var_t[dv]
            if eq_var_t !== nothing # if we already have `v_t`
                eq_idx, v_t = eq_var_t
                push!(remove_eqs, eq_idx)
                @show v_t, fullvars[dv], fullvars[v_t]
                substitute_vars!(graph, ((dv => v_t),), idx_buffer, sub_callback!; exclude = eq_idx)
            else
                v_t = diff2term(unwrap(dd))
                for eq in 𝑑neighbors(graph, dv)
                    neweqs[eq] = substitute(neweqs[eq], fullvars[dv] => v_t)
                end
                fullvars[dv] = v_t
            end
            # update the structural information
            diff_to_var[dv] = nothing
        end
    end
    @info "" fullvars

    # `SelectedState` information is no longer needed past here. State selection
    # is done. All non-differentiated variables are algebraic variables, and all
    # variables that appear differentiated are differential variables.

    ### extract partition information
    is_solvable(eq, iv) = isa(eq, Int) && BipartiteEdge(eq, iv) in solvable_graph

    solved_equations = Int[]
    solved_variables = Int[]

    # if var is like D(x)
    isdiffvar = let diff_to_var = diff_to_var
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

    var_to_idx = Dict{Any, Int}(reverse(en) for en in enumerate(fullvars))
    iv = get_iv(state.sys)
    D = Differential(iv)
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

            has_x_t = false
            x_t_idx::Union{Nothing, Int} = nothing
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
                var_eq_matching[dx_idx] = unassigned

                for eq in 𝑑neighbors(graph, dx_idx)
                    vs = 𝑠neighbors(graph, eq)
                    length(vs) == 2 || continue
                    maybe_x_t_idx = vs[1] == dx_idx ? vs[2] : vs[1]
                    maybe_x_t = fullvars[maybe_x_t_idx]
                    difference = (neweqs[eq].lhs - neweqs[eq].rhs) - (dx - maybe_x_t)
                    # if `eq` is in the form of `D(x) ~ x_t`
                    if ModelingToolkit._iszero(difference)
                        x_t_idx = maybe_x_t_idx
                        x_t = maybe_x_t
                        eq_idx = eq
                        push!(order_lowering_eqs, eq_idx)
                        has_x_t = true
                        break
                    end
                end
            end

            if x_t_idx === nothing
                x_t = ModelingToolkit.lower_varname(ogx, iv, o)
                push!(fullvars, x_t)
                x_t_idx = add_vertex!(var_to_diff)
                add_vertex!(graph, DST)
                add_vertex!(solvable_graph, DST)
                @assert x_t_idx == ndsts(graph) == length(fullvars)
                push!(var_eq_matching, unassigned)
            end
            x_t_idx::Int

            if !has_x_t
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

            end
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
                substitute_vars!(graph, subidx, idx_buffer, sub_callback!; exclude = order_lowering_eqs)
                substitute_vars!(solvable_graph, subidx, idx_buffer; exclude = order_lowering_eqs)
            end
        end
        empty!(subinfo)
        empty!(subs)
    end

    @info "" neweqs

    # Rewrite remaining equations in terms of solved variables
    function to_mass_matrix_form(ieq)
        eq = neweqs[ieq]
        if !(eq.lhs isa Number && eq.lhs == 0)
            eq = 0 ~ eq.rhs - eq.lhs
        end
        rhs = eq.rhs
        if rhs isa Symbolic
            # Check if the RHS is solvable in all state derivatives and if those
            # the linear terms for them are all zero. If so, move them to the
            # LHS.
            dterms = [var for var in 𝑠neighbors(graph, ieq) if isdiffvar(var)]
            length(dterms) == 0 && return 0 ~ rhs
            new_rhs = rhs
            new_lhs = 0
            for iv in dterms
                var = fullvars[iv]
                # 0 ~ a * D(x) + b
                # D(x) ~ -b/a
                a, b, islinear = linear_expansion(new_rhs, var)
                au = unwrap(a)
                if !islinear
                    return 0 ~ rhs
                end
                new_lhs += var
                new_rhs = -b / a
            end
            return new_lhs ~ new_rhs
        else # a number
            if abs(rhs) > 100eps(float(rhs))
                @warn "The equation $eq is not consistent. It simplifed to 0 == $rhs."
            end
            return nothing
        end
    end

    diffeq_idxs = BitSet()
    final_eqs = Equation[]
    var_rename = zeros(Int, length(var_eq_matching))
    removed_eqs = Int[]
    removed_vars = Int[]
    subeqs = Equation[]
    idx = 0
    # Solve solvable equations
    for (iv, ieq) in enumerate(var_eq_matching)
        if is_solvable(ieq, iv)
            # We don't solve differential equations, but we will need to try to
            # convert it into the mass matrix form.
            # We cannot solve the differential variable like D(x)
            if isdiffvar(iv)
                # TODO: what if `to_mass_matrix_form(ieq)` returns `nothing`?
                push!(final_eqs, to_mass_matrix_form(ieq))
                push!(diffeq_idxs, ieq)
                var_rename[iv] = (idx += 1)
                continue
            end
            eq = neweqs[ieq]
            var = fullvars[iv]
            residual = eq.lhs - eq.rhs
            a, b, islinear = linear_expansion(residual, var)
            # 0 ~ a * var + b
            # var ~ -b/a
            if ModelingToolkit._iszero(a)
                push!(removed_eqs, ieq)
                push!(removed_vars, iv)
            else
                rhs = -b/a
                neweq = var ~ simplify ? Symbolics.simplify(rhs) : rhs
                push!(subeqs, neweq)
                push!(solved_equations, ieq)
                push!(solved_variables, iv)
            end
            var_rename[iv] = -1
        else
            var_rename[iv] = (idx += 1)
        end
    end
    @info "" fullvars
    @show fullvars[solved_variables]
    @show fullvars[removed_vars]

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
                for (i, j) in enumerate(toporder)]
    end

    # TODO: BLT sorting
    # Rewrite remaining equations in terms of solved variables
    solved_eq_set = BitSet(solved_equations)
    for ieq in 1:length(neweqs)
        (ieq in diffeq_idxs || ieq in solved_eq_set) && continue
        maybe_eq = to_mass_matrix_form(ieq)
        @show maybe_eq, neweqs[ieq]
        maybe_eq === nothing || push!(final_eqs, maybe_eq)
    end
    neweqs = final_eqs

    # Contract the vertices in the structure graph to make the structure match
    # the new reality of the system we've just created.
    #
    # TODO: fix ordering and remove equations
    graph = contract_variables(graph, var_eq_matching, solved_variables)

    # Update system
    solved_variables_set = BitSet(solved_variables)
    active_vars = setdiff(BitSet(1:length(fullvars)), solved_variables_set)
    new_var_to_diff = complete(DiffGraph(length(active_vars)))
    idx = 0
    for (v, d) in enumerate(var_to_diff)
        v′ = var_rename[v]
        (v′ > 0 && d !== nothing) || continue
        d′ = var_rename[d]
        new_var_to_diff[v′] = d′ > 0 ? d′ : nothing
    end

    @set! state.structure.graph = graph
    # Note that `eq_to_diff` is not updated
    @set! state.structure.var_to_diff = new_var_to_diff
    @set! state.fullvars = [v for (i, v) in enumerate(fullvars) if i in active_vars]

    sys = state.sys
    @set! sys.eqs = neweqs
    @set! sys.states = [fullvars[i] for i in active_vars if diff_to_var[i] === nothing]
    @set! sys.observed = [observed(sys); subeqs]
    @set! sys.substitutions = Substitutions(subeqs, deps)
    @set! state.sys = sys
    @set! sys.tearing_state = state

    return invalidate_cache!(sys)
end

function tearing(state::TearingState; kwargs...)
    state.structure.solvable_graph === nothing && find_solvables!(state; kwargs...)
    complete!(state.structure)
    @unpack graph, solvable_graph = state.structure
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
