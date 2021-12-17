"""
    uneven_invmap(n::Int, list)

returns an uneven inv map with length `n`.
"""
function uneven_invmap(n::Int, list)
    rename = zero(Int, n)
    for (i, v) in enumerate(list)
        rename[v] = i
    end
    return rename
end

# N.B. assumes `slist` and `dlist` are unique
function substitution_graph(graph, slist, dlist, var_eq_matching)
    ns = length(slist)
    nd = length(dlist)
    ns == nd || error("internal error")
    newgraph = BipartiteGraph(ns, nd)
    erename = uneven_invmap(nsrc(graph), slist)
    vrename = uneven_invmap(ndst(graph), dlist)
    for e in 𝑠vertices(graph)
        ie = erename[e]
        ie == 0 && continue
        for v in 𝑠neighbors(graph, e)
            iv = vrename[v]
            iv == 0 && continue
            add_edge!(newgraph, ie, iv)
        end
    end

    newmatching = zero(slist)
    for (v, e) in enumerate(var_eq_matching)
        iv = vrename[v]
        ie = erename[e]
        iv == 0 && continue
        ie == 0 && error("internal error")
        newmatching[iv] = ie
    end

    return newgraph, newmatching
end

function tearing_sub(expr, dict, s)
    expr = ModelingToolkit.fixpoint_sub(expr, dict)
    s ? simplify(expr) : expr
end

function tearing_substitution(sys::AbstractSystem; simplify=false)
    (has_substitutions(sys) && !isnothing(get_substitutions(sys))) || return sys
    subs = get_substitutions(sys)
    neweqs = map(equations(sys)) do eq
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
                error("tearing failled because the system is singular")
            end
        end
        eq
    end
    @set! sys.eqs = neweqs
end

function solve_equation(eq, var, simplify)
    rhs = value(solve_for(eq, var; simplify=simplify, check=false))
    occursin(var, rhs) && error("solving $rhs for [$var] failed")
    var ~ rhs
end

function normalize_equation(eq)
    if !isdiffeq(eq)
        if !(eq.lhs isa Number && eq.lhs == 0)
            eq = 0 ~ eq.rhs - eq.lhs
        end
    end
    eq
end

function tearing_reassemble(sys, var_eq_matching; simplify=false)
    s = structure(sys)
    @unpack fullvars, solvable_graph, graph = s

    eqs = equations(sys)

    ### extract partition information
    is_solvable(eq, iv) = eq !== unassigned && BipartiteEdge(eq, iv) in solvable_graph

    solved_equations = Int[]
    solved_variables = Int[]

    # Solve solvable equations
    for (iv, ieq) in enumerate(var_eq_matching)
        #is_solvable(ieq, iv) || continue
        is_solvable(ieq, iv) || error("unreachable reached")
        push!(solved_equations, ieq); push!(solved_variables, iv)
    end
    subgraph, submatching = substitution_graph(graph, slist, dlist, var_eq_matching)
    toporder = topological_sort_by_dfs(DiCMOBiGraph{true}(subgraph, submatching))
    substitutions = [solve_equation(eqs[solved_equations[i]], fullvars[solved_variables[i]], simplify) for i in toporder]

    # Rewrite remaining equations in terms of solved variables

    solved_eq_set = BitSet(solved_equations)
    neweqs = Equation[normalize_equation(eqs[ieq]) for ieq in 1:length(eqs) if !(ieq in solved_eq_set)]

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
    @set! sys.substitutions = substitutions
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
