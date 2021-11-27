function tearing_sub(expr, dict, s)
    expr = ModelingToolkit.fixpoint_sub(expr, dict)
    s ? simplify(expr) : expr
end

function tearing_reassemble(sys, var_eq_matching; simplify=false)
    s = structure(sys)
    @unpack fullvars, solvable_graph, graph = s

    eqs = equations(sys)

    ### extract partition information
    function solve_equation(ieq, iv)
        var = fullvars[iv]
        eq = eqs[ieq]
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
    is_solvable(eq, iv) = eq !== unassigned && BipartiteEdge(eq, iv) in solvable_graph

    solved_equations = Int[]
    solved_variables = Int[]

    # Solve solvable equations
    for (iv, ieq) in enumerate(var_eq_matching);
        is_solvable(ieq, iv) || continue
        push!(solved_equations, ieq); push!(solved_variables, iv)
    end

    solved = Dict(solve_equation(ieq, iv) for (ieq, iv) in zip(solved_equations, solved_variables))
    obseqs = [var ~ rhs for (var, rhs) in solved]

    # Rewrite remaining equations in terms of solved variables
    function substitute_equation(ieq)
        eq = eqs[ieq]
        if isdiffeq(eq)
            return eq.lhs ~ tearing_sub(eq.rhs, solved, simplify)
        else
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
    end

    neweqs = Any[substitute_equation(ieq) for ieq in 1:length(eqs) if !(ieq in solved_equations)]
    filter!(!isnothing, neweqs)

    # Contract the vertices in the structure graph to make the structure match
    # the new reality of the system we've just created.
    graph = contract_variables(graph, var_eq_matching, solved_variables)

    # Update system
    active_vars = setdiff(BitSet(1:length(fullvars)), solved_variables)
    active_eqs = setdiff(BitSet(1:length(s.algeqs)), solved_equations)

    @set! s.graph = graph
    @set! s.fullvars = [v for (i, v) in enumerate(fullvars) if i in active_vars]
    @set! s.vartype = [v for (i, v) in enumerate(s.vartype) if i in active_vars]
    @set! s.algeqs = [e for (i, e) in enumerate(s.algeqs) if i in active_eqs]

    @set! sys.structure = s
    @set! sys.eqs = neweqs
    @set! sys.states = [s.fullvars[idx] for idx in 1:length(s.fullvars) if !isdervar(s, idx)]
    @set! sys.observed = [observed(sys); obseqs]
    return sys
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
    tear_graph_modia(graph, solvable_graph;
        varfilter=var->isalgvar(s, var), eqfilter=eq->s.algeqs[eq])
end
