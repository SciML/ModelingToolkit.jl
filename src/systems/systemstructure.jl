struct DiffGraph <: Graphs.AbstractGraph{Int}
    primal_to_diff::Vector{Union{Int, Nothing}}
    diff_to_primal::Union{Nothing, Vector{Union{Int, Nothing}}}
end
function DiffGraph(n::Integer, with_badj::Bool = false)
    DiffGraph(Union{Int, Nothing}[nothing for _ in 1:n],
        with_badj ? Union{Int, Nothing}[nothing for _ in 1:n] : nothing)
end
function Graphs.edges(dg::DiffGraph)
    (i => v for (i, v) in enumerate(dg.primal_to_diff) if v !== nothing)
end
Graphs.nv(dg::DiffGraph) = length(dg.primal_to_diff)
Base.getindex(dg::DiffGraph, var::Integer) = dg.primal_to_diff[var]
function complete(dg::DiffGraph)
    diff_to_primal = Union{Int, Nothing}[nothing for _ in 1:length(dg.primal_to_diff)]
    return DiffGraph(dg.primal_to_diff, diff_to_primal)
end
function invview(dg::DiffGraph)
end
abstract type TransformationState{T} end
abstract type AbstractTearingState{T} <: TransformationState{T} end
get_fullvars(ts::TransformationState) = ts.fullvars
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
function complete!(s::SystemStructure)
    s
end
mutable struct TearingState{T <: AbstractSystem} <: AbstractTearingState{T}
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
function is_time_dependent_parameter(p, allps, iv)
           (operation(p) === getindex &&
            (args = arguments(p); length(args)) == 1 && isequal(only(args), iv))
end
function extract_top_level_statemachines(sys::System)
    return sys, System[]
end
function TearingState(sys; check = true, sort_eqs = true)
    sys = flatten(sys)
    ivs = independent_variables(sys)
    iv = length(ivs) == 1 ? ivs[1] : nothing
    eqs = flatten_equations(equations(sys))
    original_eqs = copy(eqs)
    param_derivative_map = Dict{SymbolicT, SymbolicT}()
    no_deriv_params = Set{SymbolicT}()
    fullvars = SymbolicT[]
    var2idx = Dict{SymbolicT, Int}()
    var_types = VariableType[]
    symbolic_incidence = Vector{SymbolicT}[]
    for (i, eq) in enumerate(eqs)
        incidence = Set{SymbolicT}()
        isalgeq = true
        push!(symbolic_incidence, collect(incidence))
    end
    dervaridxs = OrderedSet{Int}()
    ndervars = length(dervaridxs)
    var_to_diff = build_var_to_diff(fullvars, ndervars, var2idx, iv)
    graph = build_incidence_graph(length(fullvars), symbolic_incidence, var2idx)
    eq_to_diff = DiffGraph(nsrcs(graph))
    return TearingState{typeof(sys)}(sys, fullvars,
        SystemStructure(complete(var_to_diff), complete(eq_to_diff),
            complete(graph), nothing, var_types, false),
        Equation[], param_derivative_map, no_deriv_params, original_eqs, Equation[], typeof(sys)[])
end
function build_var_to_diff(fullvars::Vector{SymbolicT}, ndervars::Int, var2idx::Dict{SymbolicT, Int}, @nospecialize(iv::Union{SymbolicT, Nothing}))
    var_to_diff = DiffGraph(length(fullvars), true)
    return var_to_diff
end
function build_incidence_graph(nvars::Int, symbolic_incidence::Vector{Vector{SymbolicT}}, var2idx::Dict{SymbolicT, Int})
    neqs = length(symbolic_incidence)
    graph = BipartiteGraph(neqs, nvars, Val(false))
end
function mtkcompile!(state::TearingState; simplify = false,
        check_consistency = true, fully_determined = true,
        kwargs...)
    sys, mm = ModelingToolkit.alias_elimination!(state; fully_determined, kwargs...)
    fullunknowns = [observables(sys); unknowns(sys)]
    @set! sys.observed = ModelingToolkit.topsort_equations(observed(sys), fullunknowns)
end
function _mtkcompile_worker!(state::TearingState{S}, sys::S, mm::SparseMatrixCLIL{T, Int};
                             kwargs...) where {S, T}
    return sys
end
