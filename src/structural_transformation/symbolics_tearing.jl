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

function tearing_reassemble(state::TearingState, var_eq_matching; simplify = false)
    fullvars = state.fullvars
    @unpack solvable_graph, var_to_diff, eq_to_diff, graph = state.structure

    neweqs = collect(equations(state))

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
    dummy_subs = Dict()
    diff_to_var = invview(var_to_diff)
    for var in 1:length(fullvars)
        diff_to_var[var] === nothing && continue
        if var_eq_matching[diff_to_var[var]] !== SelectedState()
            v = fullvars[var]
            # convert `D(x)` to `x_t` (don't rely on the specific spelling of
            # the name)
            dummy_subs[v] = fullvars[var] = diff2term(unwrap(v))
            # update the structural information
            diff_to_var[var] = nothing
        end
    end
    if !isempty(dummy_subs)
        neweqs = map(neweqs) do eq
            0 ~ tearing_sub(eq.rhs - eq.lhs, dummy_subs, simplify)
        end
    end

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
    # If one derivative variable `D(x)` are unsolvable in all the equations it
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
    iv = independent_variable(state.sys)
    D = Differential(iv)
    nvars = ndsts(graph)
    processed = falses(nvars)
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
        # ascend to the top of differentiation chain
        while true
            if !isempty(𝑑neighbors(graph, v))
                order = level
            end
            var_to_diff[v] === nothing && break
            processed[v] = true
            v = var_to_diff[v]
            level += 1
        end

        # `diffvar` is a order `order` variable
        order > 1 || continue

        # add `D(t) ~ x_t` etc
        subs = Dict()
        ogx = x = fullvars[diffvar] # x
        ogidx = xidx = diffvar
        for o in 1:order
            # D(x) ~ x_t
            x_t = ModelingToolkit.lower_varname(ogx, iv, o)
            dx = D(x)
            ogidx = var_to_diff[ogidx]

            x_t_idx = get(var_to_idx, x_t, nothing)
            x_t_idx !== nothing && continue

            # TODO: check x_t is legal when `x_t_idx isa Int`
            push!(fullvars, x_t)
            x_t_idx = add_vertex!(var_to_diff)
            add_vertex!(graph, DST)
            add_vertex!(solvable_graph, DST)
            @assert x_t_idx == ndsts(graph) == length(fullvars)
            push!(var_eq_matching, unassigned)

            dx_idx = get(var_to_idx, dx, nothing)
            if dx_idx === nothing
                push!(fullvars, dx)
                dx_idx = add_vertex!(var_to_diff)
                add_vertex!(graph, DST)
                add_vertex!(solvable_graph, DST)
                @assert dx_idx == ndsts(graph) == length(fullvars)
                push!(var_eq_matching, SelectedState())
            end
            add_edge!(var_to_diff, xidx, dx_idx)

            push!(neweqs, dx ~ x_t)
            eq_idx = add_vertex!(eq_to_diff)
            add_vertex!(graph, SRC)
            add_vertex!(solvable_graph, SRC)
            @assert eq_idx == nsrcs(graph) == length(neweqs)

            add_edge!(solvable_graph, eq_idx, x_t_idx)
            add_edge!(solvable_graph, eq_idx, dx_idx)
            add_edge!(graph, eq_idx, x_t_idx)
            add_edge!(graph, eq_idx, dx_idx)

            o > 1 && for eq in 𝑑neighbors(graph, ogidx)
                eq == eq_idx && continue # skip the equation that we just added
                rem_edge!(graph, eq, ogidx)
                BipartiteEdge(eq, ogidx) in solvable_graph && rem_edge!(solvable_graph, eq, ogidx)
                # TODO: what about `solvable_graph`?
                add_edge!(graph, eq, x_t_idx)
                subs[fullvars[ogidx]] = x_t
                neweqs[eq] = substitute(neweqs[eq], subs)
                empty!(subs)
            end

            # D(x_t) ~ x_tt
            x = x_t
            xidx = x_t_idx
        end
    end

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
    diffeqs = Equation[]
    # Solve solvable equations
    for (iv, ieq) in enumerate(var_eq_matching)
        is_solvable(ieq, iv) || continue
        # We don't solve differential equations, but we will need to try to
        # convert it into the mass matrix form.
        # We cannot solve the differential variable like D(x)
        if isdiffvar(iv)
            push!(diffeqs, to_mass_matrix_form(ieq))
            push!(diffeq_idxs, ieq)
            continue
        end
        push!(solved_equations, ieq)
        push!(solved_variables, iv)
    end

    if isempty(solved_equations)
        subeqs = Equation[]
        deps = Vector{Int}[]
    else
        subgraph = substitution_graph(graph, solved_equations, solved_variables,
                                      var_eq_matching)
        toporder = topological_sort_by_dfs(subgraph)
        subeqs = Equation[solve_equation(neweqs[solved_equations[i]],
                                         fullvars[solved_variables[i]],
                                         simplify) for i in toporder]
        # Find the dependency of solved variables. We will need this for ODAEProblem
        invtoporder = invperm(toporder)
        deps = [Int[invtoporder[n]
                    for n in neighborhood(subgraph, j, Inf, dir = :in) if n != j]
                for (i, j) in enumerate(toporder)]
    end

    # TODO: BLT sorting
    # Rewrite remaining equations in terms of solved variables
    solved_eq_set = BitSet(solved_equations)
    neweqs = Equation[to_mass_matrix_form(ieq)
                      for ieq in 1:length(neweqs)
                      if !(ieq in diffeq_idxs || ieq in solved_eq_set)]
    filter!(!isnothing, neweqs)
    prepend!(neweqs, diffeqs)

    # Contract the vertices in the structure graph to make the structure match
    # the new reality of the system we've just created.
    graph = contract_variables(graph, var_eq_matching, solved_variables)

    # Update system
    active_vars = setdiff(BitSet(1:length(fullvars)), solved_variables)

    @set! state.structure.graph = graph
    @set! state.fullvars = [v for (i, v) in enumerate(fullvars) if i in active_vars]

    sys = state.sys
    @set! sys.eqs = neweqs
    function isstatediff(i)
        var_eq_matching[i] !== SelectedState() && invview(var_to_diff)[i] !== nothing &&
            var_eq_matching[invview(var_to_diff)[i]] === SelectedState()
    end
    @set! sys.states = [fullvars[i] for i in active_vars if !isstatediff(i)]
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
