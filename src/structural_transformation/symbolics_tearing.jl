function tearing_sub(expr, dict, s)
    expr = ModelingToolkit.fixpoint_sub(expr, dict)
    s ? simplify(expr) : expr
end

function var_derivative!(ts::TearingState{ODESystem}, v::Int)
    sys = ts.sys
    s = ts.structure
    D = Differential(get_iv(sys))
    add_vertex!(s.solvable_graph, DST)
    push!(ts.fullvars, D(ts.fullvars[v]))
end

function eq_derivative!(ts::TearingState{ODESystem}, ieq::Int)
    sys = ts.sys
    s = ts.structure
    D = Differential(get_iv(sys))
    eq = equations(ts)[ieq]
    eq = ModelingToolkit.expand_derivatives(0 ~ D(eq.rhs - eq.lhs))
    add_vertex!(s.solvable_graph, SRC)
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
    find_eq_solvables!(ts, eq_diff; may_be_zero=true, allow_symbolic=true)
end

function tearing_reassemble(state::TearingState, var_eq_matching; simplify=false)
    fullvars = state.fullvars
    @unpack solvable_graph, var_to_diff, eq_to_diff, graph = state.structure

    neweqs = collect(equations(state))

    ### Replace derivatives of non-selected states by dumy derivatives
    dummy_subs = Dict()
    for var = 1:length(fullvars)
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
    function solve_equation(ieq, iv)
        var = fullvars[iv]
        eq = neweqs[ieq]
        rhs = value(solve_for(eq, var; simplify=simplify, check=false))

        if var in vars(rhs)
            # Usually we should be done here, but if we don't simplify we can get in
            # trouble, so try our best to still solve for rhs
            if !simplify
                rhs = SymbolicUtils.polynormalize(rhs)
            end

            # Since we know `eq` is linear wrt `var`, so the round off must be a
            # linear term. We can correct the round off error by a linear
            # correction.
            rhs -= expand_derivatives(Differential(var)(rhs))*var
            (var in vars(rhs)) && throw(EquationSolveErrors(eq, var, rhs))
        end
        var => rhs
    end
    is_solvable(eq, iv) = isa(eq, Int) && BipartiteEdge(eq, iv) in solvable_graph

    solved_equations = Int[]
    solved_variables = Int[]

    # Solve solvable equations
    for (iv, ieq) in enumerate(var_eq_matching);
        is_solvable(ieq, iv) || continue
        push!(solved_equations, ieq); push!(solved_variables, iv)
    end

    isdiffvar(var) = invview(var_to_diff)[var] !== nothing && var_eq_matching[invview(var_to_diff)[var]] === SelectedState()
    solved = Dict(solve_equation(ieq, iv) for (ieq, iv) in zip(solved_equations, solved_variables))
    obseqs = [var ~ rhs for (var, rhs) in solved]

    # Rewrite remaining equations in terms of solved variables
    function substitute_equation(ieq)
        eq = neweqs[ieq]
        if !(eq.lhs isa Number && eq.lhs == 0)
            eq = 0 ~ eq.rhs - eq.lhs
        end
        rhs = tearing_sub(eq.rhs, solved, simplify)
        if rhs isa Symbolic
            # Check if the rhs is solvable in all state derivatives and if those
            # the linear terms for them are all zero. If so, move them to the
            # LHS.
            dterms = [var for var in 𝑠neighbors(graph, ieq) if isdiffvar(var)]
            length(dterms) == 0 && return 0 ~ rhs
            new_rhs = rhs
            new_lhs = 0
            nnegative = 0
            for iv in dterms
                var = fullvars[iv]
                a, b, islinear = linear_expansion(new_rhs, var)
                au = unwrap(a)
                if !islinear || (au isa Symbolic) || isinput(var) || !(au isa Number)
                    return 0 ~ rhs
                end
                if -au < 0
                    nnegative += 1
                end
                new_lhs -= a*var
                new_rhs = b
            end
            # If most of the terms are negative, just multiply through by -1
            # to make the equations looks slightly nicer.
            if nnegative > div(length(dterms), 2)
                new_lhs = -new_lhs
                new_rhs = -new_rhs
            end
            return new_lhs ~ new_rhs
        else # a number
            if abs(rhs) > 100eps(float(rhs))
                @warn "The equation $eq is not consistent. It simplifed to 0 == $rhs."
            end
            return nothing
        end
    end

    diffeqs = [fullvars[iv] ~ tearing_sub(solved[fullvars[iv]], solved, simplify) for iv in solved_variables if isdiffvar(iv)]
    neweqs = Any[substitute_equation(ieq) for ieq in 1:length(neweqs) if !(ieq in solved_equations)]
    filter!(!isnothing, neweqs)
    prepend!(neweqs, diffeqs)

    # Update system
    active_vars = setdiff(BitSet(1:length(fullvars)), solved_variables)

    sys = state.sys
    @set! sys.eqs = neweqs
    isstatediff(i) = var_eq_matching[i] !== SelectedState() && invview(var_to_diff)[i] !== nothing && var_eq_matching[invview(var_to_diff)[i]] === SelectedState()
    @set! sys.states = [fullvars[i] for i in active_vars if !isstatediff(i)]
    @set! sys.observed = [observed(sys); obseqs]
    return sys
end

function tearing(state::TearingState)
    state.structure.solvable_graph === nothing && find_solvables!(state)
    complete!(state.structure)
    @unpack graph, solvable_graph = state.structure
    algvars = BitSet(findall(v->isalgvar(state.structure, v), 1:ndsts(graph)))
    aeqs = algeqs(state.structure)
    var_eq_matching = Matching{Union{Unassigned, SelectedState}}(tear_graph_modia(graph, solvable_graph;
        varfilter=var->var in algvars, eqfilter=eq->eq in aeqs))
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
function tearing(sys::AbstractSystem; simplify=false)
    state = TearingState(sys)
    var_eq_matching = tearing(state)
    tearing_reassemble(state, var_eq_matching; simplify=simplify)
end

"""
    tearing(sys; simplify=false)

Perform partial state selection and tearing.
"""
function partial_state_selection(sys; simplify=false)
    state = TearingState(sys)
    var_eq_matching = partial_state_selection_graph!(state)

    tearing_reassemble(state, var_eq_matching; simplify=simplify)
end
