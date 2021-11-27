function tearing_sub(expr, dict, s)
    expr = ModelingToolkit.fixpoint_sub(expr, dict)
    s ? simplify(expr) : expr
end

function tearing_reassemble(sys, var_eq_matching, eq_to_diff=nothing; simplify=false)
    s = structure(sys)
    @unpack fullvars, solvable_graph, var_to_diff, graph = s

    eqs = equations(sys)

    ### Add the differentiated equations and variables
    D = Differential(get_iv(sys))
    if length(fullvars) != length(var_to_diff)
        for i = (length(fullvars)+1):length(var_to_diff)
            push!(fullvars, D(fullvars[invview(var_to_diff)[i]]))
        end
    end

    ### Add the differentiated equations
    neweqs = copy(eqs)
    if eq_to_diff !== nothing
        eq_to_diff = complete(eq_to_diff)
        for i = (length(eqs)+1):length(eq_to_diff)
            eq = neweqs[invview(eq_to_diff)[i]]
            push!(neweqs, ModelingToolkit.expand_derivatives(0 ~ D(eq.rhs - eq.lhs)))
        end

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
            return 0 ~ rhs
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

    # Contract the vertices in the structure graph to make the structure match
    # the new reality of the system we've just created.
    graph = contract_variables(graph, var_eq_matching, solved_variables)

    # Update system
    active_vars = setdiff(BitSet(1:length(fullvars)), solved_variables)
    active_eqs = setdiff(BitSet(1:length(s.algeqs)), solved_equations)

    @set! s.graph = graph
    @set! s.fullvars = [v for (i, v) in enumerate(fullvars) if i in active_vars]
    @set! s.var_to_diff = DiffGraph(Union{Int, Nothing}[v for (i, v) in enumerate(s.var_to_diff) if i in active_vars])
    @set! s.vartype = [v for (i, v) in enumerate(s.vartype) if i in active_vars]
    @set! s.algeqs = [e for (i, e) in enumerate(s.algeqs) if i in active_eqs]

    @set! sys.structure = s
    @set! sys.eqs = neweqs
    @set! sys.states = [fullvars[i] for i in active_vars]
    @set! sys.observed = [observed(sys); obseqs]
    return sys
end

function init_for_tearing(sys)
    s = get_structure(sys)
    if !(s isa SystemStructure)
        sys = initialize_system_structure(sys)
        s = structure(sys)
    end
    find_solvables!(sys)
    @unpack graph, solvable_graph = s
    graph = complete(graph)
    @set! s.graph = graph
    @set! sys.structure = s
    return sys
end

function tear_graph(sys)
    s = structure(sys)
    @unpack graph, solvable_graph = s
    var_eq_matching = Matching{Union{Unassigned, SelectedState}}(tear_graph_modia(graph, solvable_graph;
        varfilter=var->isalgvar(s, var), eqfilter=eq->s.algeqs[eq]))
    for var in 1:ndsts(graph)
        if !isalgvar(s, var)
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
function tearing(sys; simplify=false)
    sys = init_for_tearing(sys)
    var_eq_matching = tear_graph(sys)

    tearing_reassemble(sys, var_eq_matching; simplify=simplify)
end

"""
    tearing(sys; simplify=false)

Perform partial state selection and tearing.
"""
function partial_state_selection(sys; simplify=false)
    sys = init_for_tearing(sys)
    sys, var_eq_matching, eq_to_diff = partial_state_selection_graph!(sys)

    tearing_reassemble(sys, var_eq_matching, eq_to_diff; simplify=simplify)
end
