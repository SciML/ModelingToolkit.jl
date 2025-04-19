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
    isscheduled::Bool

    function System(
            tag, eqs, noise_eqs, jumps, constraints, costs, consolidate, unknowns, ps,
            brownians, iv, observed, parameter_dependencies, var_to_name, name, description,
            defaults, guesses, systems, initialization_eqs, continuous_events, discrete_events,
            connector_type, assertions = Dict{BasicSymbolic, String}(),
            metadata = nothing, gui_metadata = nothing,
            is_dde = false, tstops = [], tearing_state = nothing, namespacing = true,
            complete = false, index_cache = nothing, ignored_connections = nothing,
            parent = nothing, isscheduled = false; checks::Union{Bool, Int} = true)
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
            parent, isscheduled)
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

    cstrs = get(kwargs, :constraints, Equation[])
    cstrunknowns, cstrps = process_constraint_system(cstrs, allunknowns, ps, iv)
    union!(allunknowns, cstrunknowns)
    union!(ps, cstrps)

    for ssys in get(kwargs, :systems, System[])
        collect_scoped_vars!(allunknowns, ps, ssys, iv)
    end

    costs = get(kwargs, :costs, nothing)
    if costs !== nothing
        costunknowns, costps = process_costs(costs, allunknowns, ps, iv)
        union!(allunknowns, costunknowns)
        union!(ps, costps)
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
Process variables in constraints of the (ODE) System.
"""
function process_constraint_system(
        constraints::Vector{Equation}, sts, ps, iv; consname = :cons)
    isempty(constraints) && return Set(), Set()

    constraintsts = OrderedSet()
    constraintps = OrderedSet()
    for cons in constraints
        collect_vars!(constraintsts, constraintps, cons, iv)
        union!(constraintsts, collect_applied_operators(cons, Differential))
    end

    # Validate the states.
    validate_vars_and_find_ps!(constraintsts, constraintps, sts, iv)

    return constraintsts, constraintps
end

"""
Process the costs for the constraint system.
"""
function process_costs(costs::Vector, sts, ps, iv)
    coststs = OrderedSet()
    costps = OrderedSet()
    for cost in costs
        collect_vars!(coststs, costps, cost, iv)
    end

    validate_vars_and_find_ps!(coststs, costps, sts, iv)
    coststs, costps
end

"""
Validate that all the variables in an auxiliary system of the (ODE) System (constraint or costs) are 
well-formed states or parameters.
 - Callable/delay variables (e.g. of the form x(0.6) should be unknowns of the system (and have one arg, etc.)
 - Callable/delay parameters should be parameters of the system

Return the set of additional parameters found in the system, e.g. in x(p) ~ 3 then p should be added as a 
parameter of the system.
"""
function validate_vars_and_find_ps!(auxvars, auxps, sysvars, iv)
    sts = sysvars

    for var in auxvars
        if !iscall(var)
            occursin(iv, var) && (var ∈ sts ||
             throw(ArgumentError("Time-dependent variable $var is not an unknown of the system.")))
        elseif length(arguments(var)) > 1
            throw(ArgumentError("Too many arguments for variable $var."))
        elseif length(arguments(var)) == 1
            if iscall(var) && operation(var) isa Differential
                var = only(arguments(var))
            end
            arg = only(arguments(var))
            operation(var)(iv) ∈ sts ||
                throw(ArgumentError("Variable $var is not a variable of the ODESystem. Called variables must be variables of the ODESystem."))

            isequal(arg, iv) || isparameter(arg) || arg isa Integer ||
                arg isa AbstractFloat ||
                throw(ArgumentError("Invalid argument specified for variable $var. The argument of the variable should be either $iv, a parameter, or a value specifying the time that the constraint holds."))

            isparameter(arg) && !isequal(arg, iv) && push!(auxps, arg)
        else
            var ∈ sts &&
                @warn "Variable $var has no argument. It will be interpreted as $var($iv), and the constraint will apply to the entire interval."
        end
    end
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

has_massactionjumps(js::System) = any(x -> x isa MassActionJump, jumps(js))
has_constantratejumps(js::System) = any(x -> x isa ConstantRateJump, jumps(js))
has_variableratejumps(js::System) = any(x -> x isa VariableRateJump, jumps(js))
# TODO: do we need this? it's kind of weird to keep
has_equations(js::System) = !isempty(equations(js))

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

function supports_initialization(sys::System)
    return isempty(jumps(sys)) && _iszero(cost(sys)) &&
           isempty(constraints(sys))
end
