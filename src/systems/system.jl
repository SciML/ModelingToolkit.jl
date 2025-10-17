struct Schedule
end
const MetadataT = Base.ImmutableDict{DataType, Any}
abstract type MutableCacheKey end
const MutableCacheT = Dict{DataType, Any}
struct System <: IntermediateDeprecationSystem
    tag::UInt
    eqs::Vector{Equation}
    noise_eqs::Union{Nothing, Vector{SymbolicT}, Matrix{SymbolicT}}
    jumps::Vector{JumpType}
    constraints::Vector{Union{Equation, Inequality}}
    costs::Vector{SymbolicT}
    consolidate::Any
    unknowns::Vector{SymbolicT}
    ps::Vector{SymbolicT}
    brownians::Vector{SymbolicT}
    iv::Union{Nothing, SymbolicT}
    observed::Vector{Equation}
    parameter_dependencies::Vector{Equation}
    var_to_name::Dict{Symbol, SymbolicT}
    name::Symbol
    description::String
    defaults::SymmapT
    guesses::SymmapT
    systems::Vector{System}
    initialization_eqs::Vector{Equation}
    continuous_events::Vector{SymbolicContinuousCallback}
    discrete_events::Vector{SymbolicDiscreteCallback}
    connector_type::Any
    assertions::Dict{SymbolicT, String}
    metadata::MetadataT
    gui_metadata::Any
    is_dde::Bool
    tstops::Vector{Any}
    inputs::OrderedSet{SymbolicT}
    outputs::OrderedSet{SymbolicT}
    tearing_state::Any
    namespacing::Bool
    complete::Bool
    index_cache::Union{Nothing, IndexCache}
    ignored_connections::Union{Nothing, Vector{Connection}}
    preface::Any
    parent::Union{Nothing, System}
    initializesystem::Union{Nothing, System}
    is_initializesystem::Bool
    is_discrete::Bool
    isscheduled::Bool
    schedule::Union{Schedule, Nothing}
    function System(
            tag, eqs, noise_eqs, jumps, constraints, costs, consolidate, unknowns, ps,
            brownians, iv, observed, parameter_dependencies, var_to_name, name, description,
            defaults, guesses, systems, initialization_eqs, continuous_events, discrete_events,
            connector_type, assertions = Dict{SymbolicT, String}(),
            metadata = MetadataT(), gui_metadata = nothing, is_dde = false, tstops = [],
            inputs = Set{SymbolicT}(), outputs = Set{SymbolicT}(),
            tearing_state = nothing, namespacing = true,
            complete = false, index_cache = nothing, ignored_connections = nothing,
            preface = nothing, parent = nothing, initializesystem = nothing,
            is_initializesystem = false, is_discrete = false, isscheduled = false,
            schedule = nothing; checks::Union{Bool, Int} = true)
        new(tag, eqs, noise_eqs, jumps, constraints, costs,
            consolidate, unknowns, ps, brownians, iv,
            observed, parameter_dependencies, var_to_name, name, description, defaults,
            guesses, systems, initialization_eqs, continuous_events, discrete_events,
            connector_type, assertions, metadata, gui_metadata, is_dde,
            tstops, inputs, outputs, tearing_state, namespacing,
            complete, index_cache, ignored_connections,
            preface, parent, initializesystem, is_initializesystem, is_discrete,
            isscheduled, schedule)
    end
end
function default_consolidate(costs, subcosts)
end
unwrap_vars(vars::AbstractArray{SymbolicT}) = vars
function unwrap_vars(vars::AbstractArray)
    for i in eachindex(vars)
    end
end
function defsdict(x::Union{AbstractDict, AbstractArray{<:Pair}})
end
function System(eqs::Vector{Equation}, iv, dvs, ps, brownians = SymbolicT[];
        constraints = Union{Equation, Inequality}[], noise_eqs = nothing, jumps = JumpType[],
        costs = SymbolicT[], consolidate = default_consolidate,
        observed = Equation[], parameter_dependencies = Equation[], defaults = SymmapT(),
        guesses = SymmapT(), systems = System[], initialization_eqs = Equation[],
        continuous_events = SymbolicContinuousCallback[], discrete_events = SymbolicDiscreteCallback[],
        connector_type = nothing, assertions = Dict{SymbolicT, String}(),
        metadata = MetadataT(), gui_metadata = nothing,
        is_dde = nothing, tstops = [], inputs = OrderedSet{SymbolicT}(),
        outputs = OrderedSet{SymbolicT}(), tearing_state = nothing,
        ignored_connections = nothing, parent = nothing,
        description = "", name = nothing, discover_from_metadata = true,
        initializesystem = nothing, is_initializesystem = false, is_discrete = false,
        preface = [], checks = true)
    if !(systems isa Vector{System})
    end
    if !isempty(parameter_dependencies)
    end
    iv = unwrap(iv)
    if iv !== nothing
    end
    if noise_eqs !== nothing
    end
    if !(inputs isa OrderedSet{SymbolicT})
    end
    if !(outputs isa OrderedSet{SymbolicT})
    end
    for subsys in systems
        for var in get_inputs(subsys)
        end
    end
    var_to_name = Dict{Symbol, SymbolicT}()
    let defaults = discover_from_metadata ? defaults : SymmapT(),
        outputs = discover_from_metadata ? outputs : OrderedSet{SymbolicT}()
    end
    continuous_events,
    discrete_events = create_symbolic_events(
        continuous_events, discrete_events)
    if is_dde === nothing
        is_dde = _check_if_dde(eqs, iv, systems)
    end
    System(Threads.atomic_add!(SYSTEM_COUNT, UInt(1)), eqs, noise_eqs, jumps, constraints,
        costs, consolidate, dvs, ps, brownians, iv, observed, Equation[],
        var_to_name, name, description, defaults, guesses, systems, initialization_eqs,
        continuous_events, discrete_events, connector_type, assertions, metadata, gui_metadata, is_dde,
        tstops, inputs, outputs, tearing_state, true, false,
        nothing, ignored_connections, preface, parent,
        initializesystem, is_initializesystem, is_discrete; checks)
end
@noinline function nonunique_subsystems(systems)
end
function System(eqs::Vector{Equation}, dvs, ps; kwargs...)
    System(eqs, nothing, dvs, ps; kwargs...)
end
function System(eqs::Vector{Equation}, iv; kwargs...)
    diffvars = OrderedSet{SymbolicT}()
    othervars = OrderedSet{SymbolicT}()
    ps = OrderedSet{SymbolicT}()
    for eq in eqs
        if !(eq.lhs isa Union{SymbolicT, Number, AbstractArray})
            if var in diffvars
            end
        end
    end
    allunknowns = union(diffvars, othervars)
    brownians = Set{SymbolicT}()
    for x in allunknowns
        if getvariabletype(x) == BROWNIAN
        end
    end
    for eq in get(kwargs, :parameter_dependencies, Equation[])
    end
    for ssys in get(kwargs, :systems, System[])
    end
    costs = get(kwargs, :costs, nothing)
    if costs !== nothing
    end
    for v in allunknowns
    end
    new_ps = gather_array_params(ps)
    noiseeqs = get(kwargs, :noise_eqs, nothing)
    if noiseeqs !== nothing
        for dv in noisedvs
        end
    end
    return System(
        eqs, iv, collect(allunknowns), collect(new_ps), collect(brownians); kwargs...)
end
function System(eqs::Vector{Equation}; kwargs...)
    return System(eqs, nothing, collect(allunknowns), collect(new_ps); kwargs...)
end
function gather_array_params(ps)
    new_ps = OrderedSet()
    for p in ps
        if iscall(p) && operation(p) === getindex
            if symbolic_has_known_size(p) && all(par[i] in ps for i in eachindex(par))
            end
        end
    end
    return new_ps
end
function process_constraint_system(
        constraints::Vector{Union{Equation, Inequality}}, sts, ps, iv; validate = true)
    isempty(constraints) && return OrderedSet{SymbolicT}(), OrderedSet{SymbolicT}()
    for cons in constraints
    end
    if validate
        if !iscall(var)
        end
    end
end
SymbolicIndexingInterface.is_time_dependent(sys::System) = get_iv(sys) !== nothing
_check_if_dde(eqs::Vector{Equation}, iv::Nothing, subsystems::Vector{System}) = false
function _check_if_dde(eqs::Vector{Equation}, iv::SymbolicT, subsystems::Vector{System})
    for eq in eqs
    end
    return false
end
function flatten(sys::System, noeqs = false)
    systems = get_systems(sys)
    isempty(systems) && return sys
    if _iszero(costs)
    end
    return System(noeqs ? Equation[] : equations(sys), get_iv(sys), unknowns(sys),
        continuous_events = continuous_events(sys),
        description = description(sys), name = nameof(sys))
end
function noise_equations_equal(sys1::System, sys2::System)
    if neqs1 === nothing && neqs2 === nothing
    end
    for eq in eqs2
    end
end
function ignored_connections_equal(sys1::System, sys2::System)
end
function NonlinearSystem(sys::System)
    if !is_time_dependent(sys)
    end
end
function SDESystem(eqs::Vector{Equation}, noise, iv; is_scalar_noise = false,
        parameter_dependencies = Equation[], kwargs...)
    if is_scalar_noise
        if !(noise isa Vector)
        end
    end
end
struct EventsInTimeIndependentSystemError <: Exception
end
function Base.showerror(io::IO, err::EventsInTimeIndependentSystemError)
    println(
        io, )
end
