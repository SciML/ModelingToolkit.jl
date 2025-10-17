struct Schedule
end
const MetadataT = Base.ImmutableDict{DataType, Any}
abstract type MutableCacheKey end
const MutableCacheT = Dict{DataType, Any}
struct System <: IntermediateDeprecationSystem
    tag::UInt
    eqs::Vector{Equation}
    constraints::Vector{Union{Equation, Inequality}}
    unknowns::Vector{SymbolicT}
    ps::Vector{SymbolicT}
    brownians::Vector{SymbolicT}
    iv::Union{Nothing, SymbolicT}
    observed::Vector{Equation}
    parameter_dependencies::Vector{Equation}
    name::Symbol
    defaults::SymmapT
    guesses::SymmapT
    systems::Vector{System}
    continuous_events::Vector{SymbolicContinuousCallback}
    discrete_events::Vector{SymbolicDiscreteCallback}
    metadata::MetadataT
    is_dde::Bool
    tstops::Vector{Any}
    inputs::OrderedSet{SymbolicT}
    outputs::OrderedSet{SymbolicT}
    tearing_state::Any
    namespacing::Bool
    complete::Bool
    index_cache::Union{Nothing, IndexCache}
    parent::Union{Nothing, System}
    initializesystem::Union{Nothing, System}
    is_initializesystem::Bool
    is_discrete::Bool
    isscheduled::Bool
    schedule::Union{Schedule, Nothing}
    function System(
            tag,
             eqs,
             constraints,
             unknowns,
             ps,
            brownians,
             iv,
             observed,
             parameter_dependencies,
             name,
            defaults,
             guesses,
             systems,
             continuous_events,
             discrete_events,
            metadata = MetadataT(),
             is_dde = false,
             tstops = [],
            inputs = Set{SymbolicT}(),
             outputs = Set{SymbolicT}(),
            tearing_state = nothing,
             namespacing = true,
            complete = false,
             index_cache = nothing,
            parent = nothing,
             initializesystem = nothing,
            is_initializesystem = false,
             is_discrete = false,
             isscheduled = false,
            schedule = nothing;
             checks::Union{Bool,Int} = true
             )
        new(
        tag,
         eqs,
         constraints,
            unknowns,
             ps,
             brownians,
             iv,
            observed,
             parameter_dependencies,
             name,
             defaults,
            guesses,
             systems,
             continuous_events,
             discrete_events,
            metadata,
             is_dde,
            tstops,
             inputs,
             outputs,
             tearing_state,
             namespacing,
            complete,
             index_cache,
            parent,
             initializesystem,
             is_initializesystem,
             is_discrete,
            isscheduled,
             schedule
             )
    end
end
unwrap_vars(vars::AbstractArray{SymbolicT}) = vars
function unwrap_vars(vars::AbstractArray)
end
function System(
    eqs::Vector{Equation},
     iv,
     dvs,
     ps,
     brownians = SymbolicT[];
        constraints = Union{Equation,
         Inequality}[],
        observed = Equation[],
         parameter_dependencies = Equation[],
         defaults = SymmapT(),
        guesses = SymmapT(),
         systems = System[],
        continuous_events = SymbolicContinuousCallback[],
         discrete_events = SymbolicDiscreteCallback[],
        metadata = MetadataT(),
        is_dde = nothing,
         tstops = [],
         inputs = OrderedSet{SymbolicT}(),
        outputs = OrderedSet{SymbolicT}(),
         tearing_state = nothing,
        parent = nothing,
        name = nothing,
         discover_from_metadata = true,
        initializesystem = nothing,
         is_initializesystem = false,
         is_discrete = false,
        checks = true
        )
    iv = unwrap(iv)
    continuous_events,
    discrete_events = create_symbolic_events(
        continuous_events, discrete_events)
    if is_dde === nothing
        is_dde = false
    end
    System(
    Threads.atomic_add!(SYSTEM_COUNT,
     UInt(1)),
     eqs,
     constraints,
        dvs,
         ps,
         brownians,
         iv,
         observed,
         Equation[],
        name,
         defaults,
         guesses,
         systems,
        continuous_events,
         discrete_events,
         metadata,
         is_dde,
        tstops,
         inputs,
         outputs,
         tearing_state,
         true,
         false,
        nothing,
         parent,
        initializesystem,
         is_initializesystem,
         is_discrete;
         checks
         )
end
function System(eqs::Vector{Equation}, dvs, ps; kwargs...)
    System(eqs, nothing, dvs, ps; kwargs...)
end
function System(eqs::Vector{Equation}, iv; kwargs...)
    diffvars = OrderedSet{SymbolicT}()
    othervars = OrderedSet{SymbolicT}()
    ps = OrderedSet{SymbolicT}()
    allunknowns = union(diffvars, othervars)
    brownians = Set{SymbolicT}()
    new_ps = gather_array_params(ps)
    return System(
        eqs, iv, collect(allunknowns), collect(new_ps), collect(brownians); kwargs...)
end
function System(eqs::Vector{Equation}; kwargs...)
    return System(eqs, nothing, collect(allunknowns), collect(new_ps); kwargs...)
end
function gather_array_params(ps)
    new_ps = OrderedSet()
    return new_ps
end
function process_constraint_system(
        constraints::Vector{Union{Equation, Inequality}}, sts, ps, iv; validate = true)
    isempty(constraints) && return OrderedSet{SymbolicT}(), OrderedSet{SymbolicT}()
end
SymbolicIndexingInterface.is_time_dependent(sys::System) = get_iv(sys) !== nothing
function flatten(sys::System, noeqs = false)
    systems = get_systems(sys)
    isempty(systems) && return sys
    return System(noeqs ? Equation[] : equations(sys), get_iv(sys), unknowns(sys),
        continuous_events = continuous_events(sys),
        name = nameof(sys))
end
