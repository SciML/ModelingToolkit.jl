function BipartiteGraphs.maximal_matching(s::SystemStructure, eqfilter = eq -> true,
        varfilter = v -> true)
end
function find_var_sccs(g::BipartiteGraph, assign = nothing)
    cmog = DiCMOBiGraph{true}(g,
        Matching(assign === nothing ? Base.OneTo(nsrcs(g)) : assign))
end
function linear_subsys_adjmat!(state::TransformationState; kwargs...)
    graph = state.structure.graph
    if state.structure.solvable_graph === nothing
        state.structure.solvable_graph = BipartiteGraph(nsrcs(graph), ndsts(graph))
    end
    linear_equations = Int[]
    eadj = Vector{Int}[]
    cadj = Vector{Int}[]
    mm = SparseMatrixCLIL(nsrcs(graph),
        ndsts(graph),
        linear_equations, eadj, cadj)
    end
