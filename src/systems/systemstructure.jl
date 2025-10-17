struct DiffGraph <: Graphs.AbstractGraph{Int}
    primal_to_diff::Vector
    diff_to_primal::Union{Nothing, Vector}
end
function foobar2 end
function foobar3 end
Base.@kwdef mutable struct SystemStructure
    var_to_diff::DiffGraph
    eq_to_diff::DiffGraph
    graph::BipartiteGraph
    solvable_graph::Union{BipartiteGraph, Nothing}
    var_types::Vector
    only_discrete::Bool
end
mutable struct TearingState{T <: AbstractSystem}
    sys::T
    fullvars::Vector
    structure::SystemStructure
    extra_eqs::Vector
    param_derivative_map::Dict
    no_deriv_params::Set
    original_eqs::Vector
    additional_observed::Vector
    statemachines::Vector
end
function extract_top_level_statemachines(sys::System)
end
function TearingState(sys; check = true, sort_eqs = true)
    sys = flatten(sys)
    eqs = equations(sys)
    var_to_diff = DiffGraph([], [])
    graph = BipartiteGraph(0, 0, Val(false))
    eq_to_diff = DiffGraph([], [])
    return TearingState(
    sys,
     SymbolicT[],
        SystemStructure(
        (var_to_diff),
         (eq_to_diff),
            graph,
             nothing,
             VariableType[],
             false
             ),
        Equation[],
         Dict(),
         Set(),
         copy(eqs),
         Equation[],
         typeof(sys)[]
         )
end
