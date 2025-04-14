struct System <: AbstractSystem
    tag::UInt
    eqs::Vector{Equation}
    # nothing - no noise
    # vector - diagonal noise
    # matrix - generic form
    # column matrix - scalar noise
    noise_eqs::Union{Nothing, AbstractVector, AbstractMatrix}
    jumps::Vector{Any}
    constraints::Vector{Union{Equation, Inequality}}
    costs::Vector{<:BasicSymbolic}
    consolidate::Any
    unknowns::Vector
    ps::Vector
    brownians::Vector
    iv::Union{Nothing, BasicSymbolic{Real}}
    observed::Vector{Equation}
    parameter_dependencies::Vector{Equation}
    var_to_name::Dict{Symbol, Any}
    name::Symbol
    description::String
    defaults::Dict
    guesses::Dict
    systems::Vector{System}
    initialization_eqs::Vector{Equation}
    continuous_events::Vector{SymbolicContinuousCallback}
    discrete_events::Vector{SymbolicDiscreteCallback}
    connector_type::Any
    assertions::Dict{BasicSymbolic, String}
    metadata::Any
    gui_metadata::Any # ?
    is_dde::Bool
    tstops::Vector{Any}
    tearing_state::Any
    namespacing::Bool
    complete::Bool
    index_cache::Union{Nothing, IndexCache}
    ignored_connections::Union{
        Nothing, Tuple{Vector{IgnoredAnalysisPoint}, Vector{IgnoredAnalysisPoint}}}
    parent::Union{Nothing, System}

    function System(
            tag, eqs, noise_eqs, jumps, constraints, costs, consolidate, unknowns, ps,
            brownians, iv, observed, parameter_dependencies, var_to_name, name, description,
            defaults, guesses, systems, initialization_eqs, continuous_events, discrete_events,
            connector_type, assertions = Dict{BasicSymbolic, String}(),
            metadata = nothing, gui_metadata = nothing,
            is_dde = false, tstops = [], tearing_state = nothing, namespacing = true,
            complete = false, index_cache = nothing, ignored_connections = nothing,
            parent = nothing; checks::Union{Bool, Int} = true)
        if (checks == true || (checks & CheckComponents) > 0) && iv !== nothing
            check_independent_variables([iv])
            check_variables(unknowns, iv)
            check_parameters(ps, iv)
            check_equations(eqs, iv)
            if noise_eqs !== nothing && size(noise_eqs, 1) != length(eqs)
                throw(IllFormedNoiseEquationsError(size(noise_eqs, 1), length(eqs)))
            end
            check_equations(equations(continuous_events), iv)
            check_subsystems(systems)
        end
        if checks == true || (checks & CheckUnits) > 0
            u = __get_unit_type(unknowns, ps, iv)
            check_units(u, eqs)
            noise_eqs !== nothing && check_units(u, noise_eqs)
            isempty(constraints) || check_units(u, constraints)
        end
        new(tag, eqs, noise_eqs, jumps, constraints, costs,
            consolidate, unknowns, ps, brownians, iv,
            observed, parameter_dependencies, var_to_name, name, description, defaults,
            guesses, systems, initialization_eqs, continuous_events, discrete_events,
            connector_type, assertions, metadata, gui_metadata, is_dde,
            tstops, tearing_state, namespacing, complete, index_cache, ignored_connections,
            parent)
    end
end

function default_consolidate(costs, subcosts)
    return sum(costs; init = 0.0) + sum(subcosts; init = 0.0)
end

function System(eqs, iv, dvs, ps, brownians = [];
        constraints = Union{Equation, Inequality}[], noise_eqs = nothing, jumps = [],
        costs = BasicSymbolic[], consolidate = default_consolidate,
        observed = Equation[], parameter_dependencies = Equation[], defaults = Dict(),
        guesses = Dict(), systems = System[], initialization_eqs = Equation[],
        continuous_events = SymbolicContinuousCallback[], discrete_events = SymbolicDiscreteCallback[],
        connector_type = nothing, assertions = Dict{BasicSymbolic, String}(),
        metadata = nothing, gui_metadata = nothing, is_dde = nothing, tstops = [],
        tearing_state = nothing, ignored_connections = nothing, parent = nothing,
        description = "", name = nothing, discover_from_metadata = true, checks = true)
    name === nothing && throw(NoNameError())

    iv = unwrap(iv)
    ps = unwrap.(ps)
    dvs = unwrap.(dvs)
    filter!(!Base.Fix2(isdelay, iv), dvs)
    brownians = unwrap.(brownians)

    if !(eqs isa AbstractArray)
        eqs = [eqs]
    end

    if noise_eqs !== nothing
        noise_eqs = unwrap.(noise_eqs)
    end

    parameter_dependencies, ps = process_parameter_dependencies(parameter_dependencies, ps)
    defaults = anydict(defaults)
    guesses = anydict(guesses)
    var_to_name = anydict()

    let defaults = discover_from_metadata ? defaults : Dict(),
        guesses = discover_from_metadata ? guesses : Dict()

        process_variables!(var_to_name, defaults, guesses, dvs)
        process_variables!(var_to_name, defaults, guesses, ps)
        process_variables!(
            var_to_name, defaults, guesses, [eq.lhs for eq in parameter_dependencies])
        process_variables!(
            var_to_name, defaults, guesses, [eq.rhs for eq in parameter_dependencies])
        process_variables!(var_to_name, defaults, guesses, [eq.lhs for eq in observed])
        process_variables!(var_to_name, defaults, guesses, [eq.rhs for eq in observed])
    end
    filter!(!(isnothing ∘ last), defaults)
    filter!(!(isnothing ∘ last), guesses)
    defaults = anydict([unwrap(k) => unwrap(v) for (k, v) in defaults])
    guesses = anydict([unwrap(k) => unwrap(v) for (k, v) in guesses])

    sysnames = nameof.(systems)
    unique_sysnames = Set(sysnames)
    if length(unique_sysnames) != length(sysnames)
        throw(NonUniqueSubsystemsError(sysnames, unique_sysnames))
    end

    continuous_events = SymbolicContinuousCallbacks(continuous_events)
    discrete_events = SymbolicDiscreteCallbacks(discrete_events)

    if iv === nothing && !isempty(continuous_events) || !isempty(discrete_events)
        throw(EventsInTimeIndependentSystemError(continuous_events, discrete_events))
    end

    if is_dde === nothing
        is_dde = _check_if_dde(eqs, iv, systems)
    end

    assertions = Dict{BasicSymbolic, String}(unwrap(k) => v for (k, v) in assertions)

    System(Threads.atomic_add!(SYSTEM_COUNT, UInt(1)), eqs, noise_eqs, jumps, constraints,
        costs, consolidate, dvs, ps, brownians, iv, observed, parameter_dependencies,
        var_to_name, name, description, defaults, guesses, systems, initialization_eqs,
        continuous_events, discrete_events, connector_type, assertions, metadata, gui_metadata, is_dde,
        tstops, tearing_state, true, false, nothing, ignored_connections, parent; checks)
end

function System(eqs, iv; kwargs...)
    iv === nothing && return System(eqs; kwargs...)
    diffvars, allunknowns, ps, eqs = process_equations(eqs, iv)
    brownians = Set()
    for x in allunknowns
        x = unwrap(x)
        if getvariabletype(x) == BROWNIAN
            push!(brownians, x)
        end
    end
    setdiff!(allunknowns, brownians)

    for eq in get(kwargs, :parameter_dependencies, Equation[])
        collect_vars!(allunknowns, ps, eq, iv)
    end

    for eq in get(kwargs, :constraints, Equation[])
        collect_vars!(allunknowns, ps, eq, iv)
    end

    for ssys in get(kwargs, :systems, System[])
        collect_scoped_vars!(allunknowns, ps, ssys, iv)
    end

    costs = get(kwargs, :costs, nothing)
    if costs !== nothing
        collect_vars!(allunknowns, ps, costs, iv)
    end

    for v in allunknowns
        isdelay(v, iv) || continue
        collect_vars!(allunknowns, ps, arguments(v)[1], iv)
    end

    new_ps = gather_array_params(ps)
    algevars = setdiff(allunknowns, diffvars)

    noiseeqs = get(kwargs, :noise_eqs, nothing)
    if noiseeqs !== nothing
        # validate noise equations
        noisedvs = OrderedSet()
        noiseps = OrderedSet()
        collect_vars!(noisedvs, noiseps, noiseeqs, iv)
        for dv in noisedvs
            dv ∈ allunknowns ||
                throw(ArgumentError("Variable $dv in noise equations is not an unknown of the system."))
        end
    end

    return System(eqs, iv, collect(Iterators.flatten((diffvars, algevars))),
        collect(new_ps), brownians; kwargs...)
end

function System(eqs; kwargs...)
    eqs = collect(eqs)

    allunknowns = OrderedSet()
    ps = OrderedSet()
    for eq in eqs
        collect_vars!(allunknowns, ps, eq, nothing)
    end
    for eq in get(kwargs, :parameter_dependencies, Equation[])
        collect_vars!(allunknowns, ps, eq, nothing)
    end
    for ssys in get(kwargs, :systems, System[])
        collect_scoped_vars!(allunknowns, ps, ssys, nothing)
    end
    costs = get(kwargs, :costs, nothing)
    if costs !== nothing
        costunknowns, costps = process_costs(costs, allunknowns, ps, nothing)
        union!(allunknowns, costunknowns)
        union!(ps, costps)
    end
    cstrs = Vector{Union{Equation, Inequality}}(get(kwargs, :constraints, []))
    for eq in cstrs
        collect_vars!(allunknowns, ps, eq, nothing)
    end

    new_ps = gather_array_params(ps)

    return System(eqs, nothing, collect(allunknowns), collect(new_ps); kwargs...)
end

function gather_array_params(ps)
    new_ps = OrderedSet()
    for p in ps
        if iscall(p) && operation(p) === getindex
            par = arguments(p)[begin]
            if Symbolics.shape(Symbolics.unwrap(par)) !== Symbolics.Unknown() &&
               all(par[i] in ps for i in eachindex(par))
                push!(new_ps, par)
            else
                push!(new_ps, p)
            end
        else
            if symbolic_type(p) == ArraySymbolic() &&
               Symbolics.shape(unwrap(p)) != Symbolics.Unknown()
                for i in eachindex(p)
                    delete!(new_ps, p[i])
                end
            end
            push!(new_ps, p)
        end
    end
    return new_ps
end

"""
    $(TYPEDSIGNATURES)

Check if a system is a (possibly implicit) discrete system. Hybrid systems are turned into
callbacks, so checking if any LHS is shifted is sufficient. If a variable is shifted in
the input equations there _will_ be a `Shift` equation in the simplified system.
"""
function is_discrete_system(sys::System)
    any(eq -> isoperator(eq.lhs, Shift), equations(sys))
end

SymbolicIndexingInterface.is_time_dependent(sys::System) = get_iv(sys) !== nothing

"""
    is_dde(sys::System)

Return a boolean indicating whether a system represents a set of delay
differential equations.
"""
is_dde(sys::System) = has_is_dde(sys) && get_is_dde(sys)

function _check_if_dde(eqs, iv, subsystems)
    is_dde = any(ModelingToolkit.is_dde, subsystems)
    if !is_dde
        vs = Set()
        for eq in eqs
            vars!(vs, eq)
            is_dde = any(vs) do sym
                isdelay(unwrap(sym), iv)
            end
            is_dde && break
        end
    end
    return is_dde
end

function flatten(sys::System, noeqs = false)
    systems = get_systems(sys)
    isempty(systems) && return sys

    return System(noeqs ? Equation[] : equations(sys), get_iv(sys), unknowns(sys),
        parameters(sys; initial_parameters = true), brownians(sys);
        jumps = jumps(sys), constraints = constraints(sys), costs = cost(sys),
        consolidate = default_consolidate, observed = observed(sys),
        parameter_dependencies = parameter_dependencies(sys), defaults = defaults(sys),
        guesses = guesses(sys), continuous_events = continuous_events(sys),
        discrete_events = discrete_events(sys), assertions = assertions(sys),
        is_dde = is_dde(sys), tstops = symbolic_tstops(sys),
        ignored_connections = ignored_connections(sys),
        initialization_eqs = initialization_equations(sys),
        # without this, any defaults/guesses obtained from metadata that were
        # later removed by the user will be re-added. Right now, we just want to
        # retain `defaults(sys)` as-is.
        discover_from_metadata = false,
        description = description(sys), name = nameof(sys))
end

"""
    $(TYPEDSIGNATURES)
"""
function check_complete(sys::System, obj)
    iscomplete(sys) || throw(SystemNotCompleteError(obj))
end

struct SystemNotCompleteError <: Exception
    obj::Any
end

function Base.showerror(io::IO, err::SystemNotCompleteError)
    print(io, """
    A completed system is required. Call `complete` or `structural_simplify` on the \
    system before creating a `$(err.obj)`.
    """)
end

struct IllFormedNoiseEquationsError <: Exception
    noise_eqs_rows::Int
    eqs_length::Int
end

function Base.showerror(io::IO, err::IllFormedNoiseEquationsError)
    print(io, """
    Noise equations are ill-formed. The number of rows much must number of drift \
    equations. `size(neqs, 1) == $(err.noise_eqs_rows) != length(eqs) == \
    $(err.eqs_length)`.
    """)
end

function NoNameError()
    ArgumentError("""
    The `name` keyword must be provided. Please consider using the `@named` macro.
    """)
end

struct NonUniqueSubsystemsError <: Exception
    names::Vector{Symbol}
    uniques::Set{Symbol}
end

function Base.showerror(io::IO, err::NonUniqueSubsystemsError)
    dupes = Set{Symbol}()
    for n in err.names
        if !(n in err.uniques)
            push!(dupes, n)
        end
        delete!(err.uniques, n)
    end
    println(io, "System names must be unique. The following system names were duplicated:")
    for n in dupes
        println(io, "  ", n)
    end
end

struct EventsInTimeIndependentSystemError <: Exception
    cevents::Vector
    devents::Vector
end

function Base.showerror(io::IO, err::EventsInTimeIndependentSystemError)
    println(io, """
    Events are not supported in time-indepent systems. Provide an independent variable to \
    make the system time-dependent or remove the events.

    The following continuous events were provided:
    $(err.cevents)

    The following discrete events were provided:
    $(err.devents) 
    """)
end
