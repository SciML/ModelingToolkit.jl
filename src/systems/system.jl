struct System <: IntermediateDeprecationSystem
    tag::UInt
    eqs::Vector{Equation}
    constraints::Vector{Union{Equation, Inequality}}
    unknowns::Vector{SymbolicT}
    ps::Vector{SymbolicT}
    iv::Union{Nothing, SymbolicT}
    observed::Vector{Equation}
    name::Symbol
    systems::Vector{System}
    continuous_events::Vector{SymbolicContinuousCallback}
    namespacing::Bool
    complete::Bool
    parent::Union{Nothing, System}
    initializesystem::Union{Nothing, System}
    is_discrete::Bool
    function System(
            tag,
             eqs,
             constraints,
             unknowns,
             ps,
             iv,
             observed,
             name,
             systems,
             continuous_events,
             namespacing = true,
            complete = false,
            parent = nothing,
             initializesystem = nothing,
             is_discrete = false;
             checks::Union{Bool,Int} = true
             )
        new(
        tag,
         eqs,
         constraints,
            unknowns,
             ps,
             iv,
            observed,
             name,
             systems,
             continuous_events,
             namespacing,
            complete,
            parent,
             initializesystem,
             is_discrete,
             )
    end
end
function System(
    eqs::Vector{Equation},
     iv,
     dvs,
     ps;
        constraints = Union{Equation,
         Inequality}[],
        observed = Equation[],
         systems = System[],
        continuous_events = SymbolicContinuousCallback[],
        name = nothing,
         is_discrete = false,
        )
    continuous_events = create_symbolic_events(
        continuous_events)
    System(
    Threads.atomic_add!(SYSTEM_COUNT,
     UInt(1)),
     eqs,
     constraints,
        dvs,
         ps,
         iv,
         observed,
        name,
         systems,
        continuous_events,
         )
end
function System(eqs::Vector{Equation}, dvs, ps; kwargs...)
    System(eqs, nothing, dvs, ps; kwargs...)
end
function System(eqs::Vector{Equation}, iv; kwargs...)
    return System(
        eqs, iv, SymbolicT[], SymbolicT[]; kwargs...)
end
SymbolicIndexingInterface.is_time_dependent(sys::System) = get_iv(sys) !== nothing
function flatten(sys::System, noeqs = false)
    systems = get_systems(sys)
    isempty(systems) && return sys
    return System(noeqs ? Equation[] : equations(sys), get_iv(sys), unknowns(sys),
        continuous_events = continuous_events(sys),
        name = nameof(sys))
end
