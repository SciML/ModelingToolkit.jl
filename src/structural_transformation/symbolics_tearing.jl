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
        e === unassigned && continue
        iv = vrename[v]
        ie = erename[e]
        iv == 0 && continue
        ie == 0 && error("internal error")
        newmatching[iv] = ie
    end

    return DiCMOBiGraph{true}(newgraph, complete(newmatching))
end

function tearing_sub(expr, dict, s)
    expr = ModelingToolkit.fixpoint_sub(expr, dict)
    s ? simplify(expr) : expr
end

function tearing_substitution(sys::AbstractSystem; simplify=false)
    empty_substitutions(sys) && return sys
    subs, = get_substitutions(sys)
    solved = Dict(eq.lhs => eq.rhs for eq in subs)
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
    @set! sys.substitutions = nothing
end

function tearing_assignments(sys::AbstractSystem)
    if empty_substitutions(sys)
        assignments = []
        deps = Int[]
        bf_states = Code.LazyState()
    else
        subs, deps = get_substitutions(sys)
        assignments = [Assignment(eq.lhs, eq.rhs) for eq in subs]
        bf_states = Code.NameState(Dict(eq.lhs => Symbol(eq.lhs) for eq in subs))
    end
    return assignments, deps, bf_states
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
        is_solvable(ieq, iv) || continue
        push!(solved_equations, ieq); push!(solved_variables, iv)
    end
    subgraph = substitution_graph(graph, solved_equations, solved_variables, var_eq_matching)
    toporder = topological_sort_by_dfs(subgraph)
    substitutions = [solve_equation(
                                    eqs[solved_equations[i]],
                                    fullvars[solved_variables[i]],
                                    simplify
                                   ) for i in toporder]
    invtoporder = invperm(toporder)
    deps = [[invtoporder[n] for n in neighborhood(subgraph, j, Inf, dir=:in) if n!=j] for (i, j) in enumerate(toporder)]

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
    @set! sys.observed = [observed(sys); substitutions]
    @set! sys.substitutions = substitutions, deps
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
