struct DiffGraph <: Graphs.AbstractGraph{Int}
    primal_to_diff::Vector{Union{Int, Nothing}}
    diff_to_primal::Union{Nothing, Vector{Union{Int, Nothing}}}
end
function DiffGraph(n::Integer, with_badj::Bool = false)
    DiffGraph(Union{Int, Nothing}[nothing for _ in 1:n],
        with_badj ? Union{Int, Nothing}[nothing for _ in 1:n] : nothing)
end
function complete(dg::DiffGraph)
    diff_to_primal = Union{Int, Nothing}[nothing for _ in 1:length(dg.primal_to_diff)]
    return DiffGraph(dg.primal_to_diff, diff_to_primal)
end
Base.@kwdef mutable struct SystemStructure
    var_to_diff::DiffGraph
    eq_to_diff::DiffGraph
    graph::BipartiteGraph{Int, Nothing}
    solvable_graph::Union{BipartiteGraph{Int, Nothing}, Nothing}
    var_types::Vector{VariableType}
    only_discrete::Bool
end
function Base.copy(structure::SystemStructure)
    SystemStructure(copy(structure.var_to_diff), copy(structure.eq_to_diff),
        var_types, structure.only_discrete)
end
mutable struct TearingState{T <: AbstractSystem}
    sys::T
    fullvars::Vector{SymbolicT}
    structure::SystemStructure
    extra_eqs::Vector{Equation}
    param_derivative_map::Dict{SymbolicT, SymbolicT}
    no_deriv_params::Set{SymbolicT}
    original_eqs::Vector{Equation}
    additional_observed::Vector{Equation}
    statemachines::Vector{T}
end
struct EquationsView{T} <: AbstractVector{Equation}
    ts::TearingState{T}
end
equations(ts::TearingState) = EquationsView(ts)
Base.size(ev::EquationsView) = (length(equations(ev.ts.sys)) + length(ev.ts.extra_eqs),)
function Base.getindex(ev::EquationsView, i::Integer)
    eqs = equations(ev.ts.sys)
    return eqs[i]
end
function extract_top_level_statemachines(sys::System)
    return sys, System[]
end
function TearingState(sys; check = true, sort_eqs = true)
    sys = flatten(sys)
    ivs = independent_variables(sys)
    iv = length(ivs) == 1 ? ivs[1] : nothing
    eqs = equations(sys)
    var_to_diff = DiffGraph(0, true)
    graph = BipartiteGraph(0, 0, Val(false))
    eq_to_diff = DiffGraph(0)
    return TearingState{typeof(sys)}(
    sys,
     SymbolicT[],
        SystemStructure(
        complete(var_to_diff),
         complete(eq_to_diff),
            complete(graph),
             nothing,
             VariableType[],
             false
             ),
        Equation[],
         Dict{SymbolicT, SymbolicT}(),
         Set{SymbolicT}(),
         copy(eqs),
         Equation[],
         typeof(sys)[]
         )
end
