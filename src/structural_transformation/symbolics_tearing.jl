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
    # remove those that doen't actually occur.
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

    ### Replace derivatives of non-selected states by dumy derivatives
    dummy_subs = Dict()
    for var in 1:length(fullvars)
        invview(var_to_diff)[var] === nothing && continue
        if var_eq_matching[invview(var_to_diff)[var]] !== SelectedState()
            fullvar = fullvars[var]
            subst_fullvar = tearing_sub(fullvar, dummy_subs, simplify)
            dummy_subs[fullvar] = fullvars[var] = diff2term(unwrap(subst_fullvar))
            var_to_diff[invview(var_to_diff)[var]] = nothing
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
    function isdiffvar(var)
        invview(var_to_diff)[var] !== nothing &&
            var_eq_matching[invview(var_to_diff)[var]] === SelectedState()
    end

    # Rewrite remaining equations in terms of solved variables
    function to_mass_matrix_form(ieq)
        eq = neweqs[ieq]
        if !(eq.lhs isa Number && eq.lhs == 0)
            eq = 0 ~ eq.rhs - eq.lhs
        end
        rhs = eq.rhs
        if rhs isa Symbolic
            # Check if the rhs is solvable in all state derivatives and if those
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
        # find the dependency of solved variables. we will need this for ODAEProblem
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

Perform index reduction and use the dummy derivative techinque to ensure that
the system is balanced.
"""
function dummy_derivative(sys, state = TearingState(sys); kwargs...)
    function jac(eqs, vars)
        symeqs = EquationsView(state)[eqs]
        Symbolics.jacobian((x -> x.rhs).(symeqs), state.fullvars[vars])
    end
    var_eq_matching = dummy_derivative_graph!(state, jac; kwargs...)
    tearing_reassemble(state, var_eq_matching)
end
