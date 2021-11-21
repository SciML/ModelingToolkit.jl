"""
    tear_graph(sys) -> sys

Tear the bipartite graph in a system. End users are encouraged to call [`structural_simplify`](@ref)
instead, which calls this function internally.
"""
function tear_graph(sys)
    find_solvables!(sys)
    s = structure(sys)
    @unpack graph, solvable_graph, var_eq_matching, scc = s

    @set! sys.structure.partitions = map(scc) do c
        ieqs = filter(eq->isalgeq(s, eq), c)
        vars = Int[var for var in invview(var_eq_matching)[ieqs] if var !== unassigned]

        ict = IncrementalCycleTracker(DiCMOBiGraph{true}(graph); dir=:in)
        SystemPartition(tearEquations!(ict, solvable_graph.fadjlist, ieqs, vars)...)
    end
    return sys
end

"""
    algebraic_equations_scc(sys)

Find strongly connected components of algebraic equations in a system.
"""
function algebraic_equations_scc(sys)
    s = get_structure(sys)
    if !(s isa SystemStructure)
        sys = initialize_system_structure(sys)
        s = structure(sys)
    end

    # skip over differential equations
    algvars = isalgvar.(Ref(s), 1:ndsts(s.graph))
    var_eq_matching = complete(matching(s, algvars, s.algeqs))
    components = find_scc(s.graph, var_eq_matching)

    @set! sys.structure.var_eq_matching = var_eq_matching
    @set! sys.structure.scc = components
    return sys
end
