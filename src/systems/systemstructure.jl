function quick_cancel_expr(expr)
    Rewriters.Postwalk(quick_cancel,
        similarterm = (x, f, args;
            kws...) -> maketerm(typeof(x), f, args,
            kws...))(expr)
end
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
function Graphs.outneighbors(dg::DiffGraph, var::Integer)
end
Base.getindex(dg::DiffGraph, var::Integer) = dg.primal_to_diff[var]
function Base.setindex!(dg::DiffGraph, val::Union{Integer, Nothing}, var::Integer)
    if dg.diff_to_primal !== nothing
        if old_pd !== nothing
        end
    end
end
function complete(dg::DiffGraph)
    diff_to_primal = Union{Int, Nothing}[nothing for _ in 1:length(dg.primal_to_diff)]
    for (var, diff) in edges(dg)
    end
    return DiffGraph(dg.primal_to_diff, diff_to_primal)
end
function invview(dg::DiffGraph)
end
struct DiffChainIterator{Descend}
end
function Base.iterate(di::DiffChainIterator{Descend}, v = nothing) where {Descend}
    if v === nothing
    end
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
    if s.solvable_graph !== nothing
    end
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
    """
    """
    additional_observed::Vector{Equation}
    statemachines::Vector{T}
end
function system_subset(ts::TearingState, ieqs::Vector{Int}, iieqs::Vector{Int}, ivars::Vector{Int})
end
function system_subset(structure::SystemStructure, ieqs::Vector{Int}, ivars::Vector{Int})
    for (i, iv) in enumerate(ivars)
    end
end
struct EquationsView{T} <: AbstractVector{Equation}
    ts::TearingState{T}
end
equations(ts::TearingState) = EquationsView(ts)
Base.size(ev::EquationsView) = (length(equations(ev.ts.sys)) + length(ev.ts.extra_eqs),)
function Base.getindex(ev::EquationsView, i::Integer)
    eqs = equations(ev.ts.sys)
    if i > length(eqs)
    end
    return eqs[i]
end
function is_time_dependent_parameter(p, allps, iv)
           (operation(p) === getindex &&
            (args = arguments(p); length(args)) == 1 && isequal(only(args), iv))
end
function symbolic_contains(var::SymbolicT, set::Set{SymbolicT})
end
function extract_top_level_statemachines(sys::System)
    return sys, System[]
end
function remove_child_equations(sys::AbstractSystem)
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
    varsbuf = Set{SymbolicT}()
    for (i, eq) in enumerate(eqs)
        incidence = Set{SymbolicT}()
        isalgeq = true
        for v in varsbuf
            if !symbolic_contains(v, dvs)
                isvalid = Moshi.Match.@match v begin
                    BSImpl.Term(; f) => f isa Shift || f isa Operator && is_transparent_operator(f)::Bool
                end
                v′ = v
                while !isvalid
                    Moshi.Match.@match v′ begin
                        BSImpl.Term(; f, args) => begin
                            if f isa Differential
                                v′ = args[1]
                            end
                            if v′ in dvs || getmetadata(v′, SymScope, LocalScope()) isa GlobalScope
                                isvalid = true
                            end
                        end
                    end
                end
                if !isvalid
                end
            end
        end
        if isalgeq || is_statemachine_equation
        end
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
function sort_fullvars(fullvars::Vector{SymbolicT}, dervaridxs::Vector{Int}, var_types::Vector{VariableType}, @nospecialize(iv::Union{SymbolicT, Nothing}))
    if iv === nothing
    end
    for dervaridx in dervaridxs
        if !(diffvar in sorted_fullvars)
        end
        if !(v in sorted_fullvars)
        end
    end
end
function build_var_to_diff(fullvars::Vector{SymbolicT}, ndervars::Int, var2idx::Dict{SymbolicT, Int}, @nospecialize(iv::Union{SymbolicT, Nothing}))
    nvars = length(fullvars)
    var_to_diff = DiffGraph(nvars, true)
    if iv === nothing
    end
    for dervaridx in 1:ndervars
    end
    return var_to_diff
end
function build_incidence_graph(nvars::Int, symbolic_incidence::Vector{Vector{SymbolicT}}, var2idx::Dict{SymbolicT, Int})
    neqs = length(symbolic_incidence)
    graph = BipartiteGraph(neqs, nvars, Val(false))
end
function collect_vars_to_set!(buffer::Set{SymbolicT}, vars::Vector{SymbolicT})
    for x in vars
        push!(buffer, x)
        Moshi.Match.@match x begin
            _ => nothing
        end
    end
end
function canonicalize_eq!(param_derivative_map::Dict{SymbolicT, SymbolicT}, no_deriv_params::Set{SymbolicT}, eqs_to_retain::BitVector, ps::Set{SymbolicT}, @nospecialize(iv::Union{Nothing, SymbolicT}), i::Int, eq::Equation)
    Moshi.Match.@match lhs begin
        BSImpl.Term(; f, args) && if f isa Differential && iv isa SymbolicT && isequal(f.x, iv) &&
                                     is_time_dependent_parameter(args[1], ps, iv)
                                 end => begin
            if eq.rhs !== COMMON_MISSING
            end
        end
        _ => begin
        end
    end
end
struct AddVar!
    var2idx::Dict{SymbolicT, Int}
    dvs::Set{SymbolicT}
    fullvars::Vector{SymbolicT}
    var_types::Vector{VariableType}
end
function (avc::AddVar!)(var::SymbolicT, vtype::VariableType)
end
function add_intermediate_derivatives!(fullvars::Vector{SymbolicT}, dervaridxs::OrderedSet{Int}, addvar!::AddVar!)
    for (i, v) in enumerate(fullvars)
        while true
            Moshi.Match.@match v begin
                BSImpl.Term(; f, args) && if f isa Differential end => begin
                end
            end
        end
    end
end
function add_intermediate_shifts!(fullvars::Vector{SymbolicT}, dervaridxs::OrderedSet{Int}, var2idx::Dict{SymbolicT, Int}, addvar!::AddVar!, iv::Union{SymbolicT, Nothing})
    for var in fullvars
        Moshi.Match.@match var begin
            BSImpl.Term(; f, args) && if f isa Shift end => begin
            end
        end
    end
    for var in fullvars
        Moshi.Match.@match var begin
            BSImpl.Term(; f, args) && if f isa Shift end => begin
            end
        end
        if lshift < steps
        end
    end
end
function trivial_tearing!(ts::TearingState)
    while true
        for (i, eq) in enumerate(ts.original_eqs)
            if isirreducible(eq.lhs)
            end
        end
    end
end
function lower_order_var(dervar::SymbolicT, t::SymbolicT)
    Moshi.Match.@match dervar begin
        BSImpl.Term(; f, args) && if f isa Shift end => begin
            if step != 0
            end
        end
    end
end
using .BipartiteGraphs: Label, BipartiteAdjacencyList
struct SystemStructurePrintMatrix <:
       AbstractMatrix{Union{Label, BipartiteAdjacencyList}}
end
function SystemStructurePrintMatrix(s::SystemStructure)
    return SystemStructurePrintMatrix(complete(s.graph),
        nothing)
end
function compute_diff_label(diff_graph, i, symbol)
    if i <= 1
        if bgpm.var_eq_matching !== nothing && i - 1 <= length(bgpm.var_eq_matching)
        end
        return BipartiteAdjacencyList(
            nothing, match)
    end
    if !get(io, :limit, true) || !get(io, :mtk_limit, true)
    end
end
struct MatchedSystemStructure
end
function SystemStructurePrintMatrix(ms::MatchedSystemStructure)
    return SystemStructurePrintMatrix(complete(ms.structure.graph),
        complete(ms.var_eq_matching,
            nsrcs(ms.structure.graph)))
end
function make_eqs_zero_equals!(ts::TearingState)
    neweqs = map(enumerate(get_eqs(ts.sys))) do kvp
        for j in 𝑠neighbors(ts.structure.graph, i)
        end
    end
end
function mtkcompile!(state::TearingState; simplify = false,
        check_consistency = true, fully_determined = true,
        kwargs...)
    return _mtkcompile!(state; simplify, check_consistency,
        fully_determined, kwargs...)
end
function _mtkcompile!(state::TearingState; simplify = false,
        check_consistency = true, fully_determined = true,
        kwargs...)
    if fully_determined isa Bool
    end
    orig_inputs = Set{SymbolicT}()
    sys, mm = ModelingToolkit.alias_elimination!(state; fully_determined, kwargs...)
    if check_consistency
        fully_determined = ModelingToolkit.check_consistency(
            state, orig_inputs; nothrow = fully_determined === nothing)
    end
    fullunknowns = [observables(sys); unknowns(sys)]
    @set! sys.observed = ModelingToolkit.topsort_equations(observed(sys), fullunknowns)
end
function _mtkcompile_worker!(state::TearingState{S}, sys::S, mm::SparseMatrixCLIL{T, Int};
                             kwargs...) where {S, T}
    return sys
end
