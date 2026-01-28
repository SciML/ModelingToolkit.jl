struct Schedule
    var_sccs::Vector{Vector{Int}}
    """
    Mapping of `Differential`s of variables to corresponding derivative expressions.
    """
    dummy_sub::Dict{SymbolicT, SymbolicT}
end

function Base.copy(sched::Schedule)
    return Schedule(copy(sched.var_sccs), copy(sched.dummy_sub))
end

const MetadataT = Base.ImmutableDict{DataType, Any}

abstract type MutableCacheKey end

const MutableCacheT = Dict{DataType, Any}

"""
    $TYPEDEF

Utility metadata key for adding miscellaneous/one-off metadata to systems.
"""
abstract type MiscSystemData end

"""
    $(TYPEDEF)

A symbolic representation of a numerical system to be solved. This is a recursive
tree-like data structure - each system can contain additional subsystems. As such,
it implements the `AbstractTrees.jl` interface to enable exploring the hierarchical
structure.

# Fields

$(TYPEDFIELDS)
"""
struct System <: IntermediateDeprecationSystem
    """
    $INTERNAL_FIELD_WARNING
    A unique integer tag for the system.
    """
    tag::UInt
    """
    The equations of the system.
    """
    eqs::Vector{Equation}
    # nothing - no noise
    # vector - diagonal noise
    # matrix - generic form
    # column matrix - scalar noise
    """
    The noise terms for each equation of the system. This field is only used for flattened
    systems. To represent noise in a hierarchical system, use brownians. In a system with
    `N` equations and `K` independent brownian variables, this should be an `N x K`
    matrix. In the special case where `N == K` and each equation has independent noise,
    this noise matrix is diagonal. Diagonal noise can be specified by providing an `N`
    length vector. If this field is `nothing`, the system does not have noise.
    """
    noise_eqs::Union{Nothing, Vector{SymbolicT}, Matrix{SymbolicT}}
    """
    Jumps associated with the system. Each jump can be a `VariableRateJump`,
    `ConstantRateJump` or `MassActionJump`. See `JumpProcesses.jl` for more information.
    """
    jumps::Vector{JumpType}
    """
    The constraints of the system. This can be used to represent the constraints in an
    optimal-control problem or boundary-value differential equation, or the constraints
    in a constrained optimization.
    """
    constraints::Vector{Union{Equation, Inequality}}
    """
    The costs of the system. This can be the cost in an optimal-control problem, or the
    loss of an optimization problem. Scalar loss values must also be provided as a single-
    element vector.
    """
    costs::Vector{SymbolicT}
    """
    A function which combines costs into a scalar value. This should take two arguments,
    the `costs` of this system and the consolidated costs of all subsystems in the order
    they are present in the `systems` field. It should return a scalar cost that combines
    all of the individual values. This defaults to a function that simply sums all cost
    values.
    """
    consolidate::Any
    """
    The variables being solved for by this system. For example, in a differential equation
    system, this contains the dependent variables.
    """
    unknowns::Vector{SymbolicT}
    """
    The parameters of the system. Parameters can either be variables that parameterize the
    problem being solved for (e.g. the spring constant of a mass-spring system) or
    additional unknowns not part of the main dynamics of the system (e.g. discrete/clocked
    variables in a hybrid ODE).
    """
    ps::Vector{SymbolicT}
    """
    The brownian variables of the system, created via `@brownians`. Each brownian variable
    represents an independent noise. A system with brownians cannot be simulated directly.
    It needs to be compiled using `mtkcompile` into `noise_eqs`.
    """
    brownians::Vector{SymbolicT}
    """
    The independent variable for a time-dependent system, or `nothing` for a time-independent
    system.
    """
    iv::Union{Nothing, SymbolicT}
    """
    Equations that compute variables of a system that have been eliminated from the set of
    unknowns by `mtkcompile`. More generally, this contains all variables that can be
    computed from the unknowns and parameters and do not need to be solved for. Such
    variables are termed as "observables". Each equation must be of the form
    `observable ~ expression` and observables cannot appear on the LHS of multiple
    equations. Equations must be sorted such that every observable appears on
    the left hand side of an equation before it appears on the right hand side of any other
    equation.
    """
    observed::Vector{Equation}
    """
    $INTERNAL_FIELD_WARNING
    A mapping from the name of a variable to the actual symbolic variable in the system.
    This is used to enable `getproperty` syntax to access variables of a system.
    """
    var_to_name::Dict{Symbol, SymbolicT}
    """
    The name of the system.
    """
    name::Symbol
    """
    An optional description for the system.
    """
    description::String
    """
    Binding relations for variables/parameters. The bound variable (key) is completely
    determined by the binding (value). Providing an initial condition for a bound variable
    is an error. Bindings for variables (ones created via `@variables` and `@discretes`)
    are treated as initial conditions.
    """
    bindings::ROSymmapT
    """
    Initial conditions for variables (unknowns/observables/parameters) which can be
    changed/overridden. When constructing a numerical problem from the system.
    """
    initial_conditions::SymmapT
    """
    Guess values for variables of a system that are solved for during initialization.
    """
    guesses::SymmapT
    """
    A list of subsystems of this system. Used for hierarchically building models.
    """
    systems::Vector{System}
    """
    Equations that must be satisfied during initialization of the numerical problem created
    from this system. For time-dependent systems, these equations are not valid after the
    initial time.
    """
    initialization_eqs::Vector{Equation}
    """
    Symbolic representation of continuous events in a dynamical system. See
    [`SymbolicContinuousCallback`](@ref).
    """
    continuous_events::Vector{SymbolicContinuousCallback}
    """
    Symbolic representation of discrete events in a dynamica system. See
    [`SymbolicDiscreteCallback`](@ref).
    """
    discrete_events::Vector{SymbolicDiscreteCallback}
    """
    $INTERNAL_FIELD_WARNING
    If this system is a connector, the type of connector it is.
    """
    connector_type::Any
    """
    A map from expressions that must be through throughout the solution process to an
    associated error message. By default these assertions cause the generated code to
    output `NaN`s if violated, but can be made to error using `debug_system`.
    """
    assertions::Dict{SymbolicT, String}
    """
    The metadata associated with this system, as a `Base.ImmutableDict`. This follows
    the same interface as SymbolicUtils.jl. Metadata can be queried and updated using
    `SymbolicUtils.getmetadata` and `SymbolicUtils.setmetadata` respectively.
    """
    metadata::MetadataT
    """
    $INTERNAL_FIELD_WARNING
    Metadata added by the `@mtkmodel` macro.
    """
    gui_metadata::Any # ?
    """
    Whether the system contains delay terms. This is inferred from the equations, but
    can also be provided explicitly.
    """
    is_dde::Bool
    """
    Extra time points for the integrator to stop at. These can be numeric values,
    or expressions of parameters and time.
    """
    tstops::Vector{Any}
    """
    $INTERNAL_FIELD_WARNING
    The list of input variables of the system.
    """
    inputs::OrderedSet{SymbolicT}
    """
    $INTERNAL_FIELD_WARNING
    The list of output variables of the system.
    """
    outputs::OrderedSet{SymbolicT}
    """
    The `TearingState` of the system post-simplification with `mtkcompile`.
    """
    tearing_state::Any
    """
    Whether the system namespaces variables accessed via `getproperty`. `complete`d systems
    do not namespace, but this flag can be toggled independently of `complete` using
    `toggle_namespacing`.
    """
    namespacing::Bool
    """
    Whether the system is marked as "complete". Completed systems cannot be used as
    subsystems.
    """
    complete::Bool
    """
    $INTERNAL_FIELD_WARNING
    For systems simplified or completed with `split = true` (the default) this contains an
    `IndexCache` which aids in symbolic indexing. If this field is `nothing`, the system is
    either not completed, or completed with `split = false`.
    """
    index_cache::Union{Nothing, IndexCache}
    """
    $INTERNAL_FIELD_WARNING
    Contains the dependency graph of bound parameters to avoid excessive duplicated work
    during code generation.
    """
    parameter_bindings_graph::Union{Nothing, ParameterBindingsGraph}
    """
    $INTERNAL_FIELD_WARNING
    Connections that should be ignored because they were removed by an analysis point
    transformation. The first element of the tuple contains all such "standard" connections
    (ones between connector systems) and the second contains all such causal variable
    connections.
    """
    ignored_connections::Union{Nothing, Vector{Connection}}
    """
    `SymbolicUtils.Code.Assignment`s to prepend to all code generated from this system.
    """
    preface::Any
    """
    After simplification with `mtkcompile`, this field contains the unsimplified system
    with the hierarchical structure. There may be multiple levels of `parent`s. The root
    parent is used for accessing variables via `getproperty` syntax.
    """
    parent::Union{Nothing, System}
    """
    A custom initialization system to use if no initial conditions are provided for the
    unknowns or observables of this system.
    """
    initializesystem::Union{Nothing, System}
    """
    Whether the current system is an initialization system.
    """
    is_initializesystem::Bool
    is_discrete::Bool
    """
    State priorities for variables. Used in structural simplification algorithms.
    """
    state_priorities::AtomicMapT{Int}
    """
    Variables marked as irreducible for simplification.
    """
    irreducibles::AtomicSetT
    """
    $INTERNAL_FIELD_WARNING
    Whether the system has been simplified by `mtkcompile`.
    """
    isscheduled::Bool
    """
    $INTERNAL_FIELD_WARNING
    The `Schedule` containing additional information about the simplified system.
    """
    schedule::Union{Schedule, Nothing}

    function System(
            tag, eqs, noise_eqs, jumps, constraints, costs, consolidate, unknowns, ps,
            brownians, iv, observed, var_to_name, name, description, bindings,
            initial_conditions, guesses, systems, initialization_eqs, continuous_events,
            discrete_events, connector_type, assertions = Dict{SymbolicT, String}(),
            metadata = MetadataT(), gui_metadata = nothing, is_dde = false, tstops = [],
            inputs = Set{SymbolicT}(), outputs = Set{SymbolicT}(),
            tearing_state = nothing, namespacing = true,
            complete = false, index_cache = nothing, parameter_bindings_graph = nothing,
            ignored_connections = nothing,
            preface = nothing, parent = nothing, initializesystem = nothing,
            is_initializesystem = false, is_discrete = false, state_priorities = AtomicMapT{Int}(),
            irreducibles = AtomicSetT(), isscheduled = false,
            schedule = nothing; checks::Union{Bool, Int} = true
        )
        if is_initializesystem && iv !== nothing
            throw(
                ArgumentError(
                    """
                    Expected initialization system to be time-independent. Found independent
                    variable $iv.
                    """
                )
            )
        end
        @assert iv === nothing || symtype(iv) === Real
        if (checks isa Bool && checks === true || checks isa Int && (checks & CheckComponents) > 0) && iv !== nothing
            check_independent_variables((iv,))
            check_variables(unknowns, iv)
            check_parameters(ps, iv)
            check_equations(eqs, iv)
            Neq = length(eqs)
            if noise_eqs isa Matrix{SymbolicT}
                N1 = size(noise_eqs, 1)
            elseif noise_eqs isa Vector{SymbolicT}
                N1 = length(noise_eqs)
            elseif noise_eqs === nothing
                N1 = Neq
            else
                error()
            end
            N1 == Neq || throw(IllFormedNoiseEquationsError(N1, Neq))
            check_equations(equations(continuous_events), iv)
            check_subsystems(systems)
        end
        if checks == true || (checks & CheckUnits) > 0
            u = __get_unit_type(unknowns, ps, iv)
            if noise_eqs === nothing
                check_units(u, eqs)
            else
                check_units(u, eqs, noise_eqs)
            end
            if iv !== nothing
                check_units(u, jumps, iv)
            end
            isempty(constraints) || check_units(u, constraints)
        end
        return new(
            tag, eqs, noise_eqs, jumps, constraints, costs,
            consolidate, unknowns, ps, brownians, iv,
            observed, var_to_name, name, description, bindings, initial_conditions,
            guesses, systems, initialization_eqs, continuous_events, discrete_events,
            connector_type, assertions, metadata, gui_metadata, is_dde,
            tstops, inputs, outputs, tearing_state, namespacing,
            complete, index_cache, parameter_bindings_graph, ignored_connections,
            preface, parent, initializesystem, is_initializesystem, is_discrete,
            state_priorities, irreducibles,
            isscheduled, schedule
        )
    end
end

_sum_costs(costs::Vector{SymbolicT}) = SU.add_worker(VartypeT, costs)
_sum_costs(costs::Vector{Num}) = SU.add_worker(VartypeT, costs)
# `reduce` instead of `sum` because the rrule for `sum` doesn't
# handle the `init` kwarg.
_sum_costs(costs::Vector) = reduce(+, costs; init = 0.0)

function default_consolidate(costs, subcosts)
    return _sum_costs(costs) + _sum_costs(subcosts)
end

unwrap_vars(x) = unwrap_vars(collect(x))
unwrap_vars(vars::AbstractArray{SymbolicT}) = vars
function unwrap_vars(vars::AbstractArray)
    result = similar(vars, SymbolicT)
    for i in eachindex(vars)
        result[i] = SU.Const{VartypeT}(vars[i])
    end
    return result
end

defsdict(x::SymmapT) = x
function defsdict(x::Union{AbstractDict, AbstractArray{<:Pair}})
    result = SymmapT()
    for (k, v) in x
        result[unwrap(k)] = SU.Const{VartypeT}(v)
    end
    return result
end

as_atomicmap(::Type{T}, x::AtomicMapT{T}) where {T} = x
function as_atomicmap(::Type{T}, x::Union{AbstractDict, AbstractArray{<:Pair}}) where {T}
    result = AtomicMapT{T}()
    for (k, v) in x
        result[unwrap(k)] = convert(T, v)
    end
    return result
end

parse_atomicset(x::AtomicSetT) = x
function parse_atomicset(vars)
    result = AtomicSetT()
    for var in vars
        push!(result, unwrap(var))
    end
    return result
end

function __get_new_tag()
    return Threads.atomic_add!(SYSTEM_COUNT, UInt(1))
end

"""
    $(TYPEDSIGNATURES)

Construct a system using the given equations `eqs`, independent variable `iv` (`nothing`)
for time-independent systems, unknowns `dvs`, parameters `ps` and brownian variables
`brownians`.

## Keyword Arguments

- `discover_from_metadata`: Whether to parse metadata of unknowns and parameters of the
  system to obtain bindings, initial conditions and/or guesses.
- `checks`: Whether to perform sanity checks on the passed values.

All other keyword arguments are named identically to the corresponding fields in
[`System`](@ref).
"""
function System(
        eqs::Vector{Equation}, iv, dvs, ps, brownians = SymbolicT[];
        constraints = Union{Equation, Inequality}[], noise_eqs = nothing, jumps = JumpType[],
        costs = SymbolicT[], consolidate = default_consolidate,
        # `@nospecialize` is only supported on the first 32 arguments. Keep this early.
        @nospecialize(preface = nothing), @nospecialize(tstops = []),
        observed = Equation[], bindings = SymmapT(), initial_conditions = SymmapT(),
        guesses = SymmapT(), systems = System[], initialization_eqs = Equation[],
        continuous_events = SymbolicContinuousCallback[], discrete_events = SymbolicDiscreteCallback[],
        connector_type = nothing, assertions = Dict{SymbolicT, String}(),
        metadata = MetadataT(), gui_metadata = nothing,
        is_dde = nothing, inputs = OrderedSet{SymbolicT}(),
        outputs = OrderedSet{SymbolicT}(), tearing_state = nothing,
        ignored_connections = nothing, parent = nothing, state_priorities = AtomicMapT{Int}(),
        irreducibles = AtomicSetT(),
        description = "", name = nothing, discover_from_metadata = true,
        initializesystem = nothing, is_initializesystem = false, is_discrete = false,
        checks = true, __legacy_defaults__ = nothing
    )
    name === nothing && throw(NoNameError())

    if __legacy_defaults__ !== nothing
        Base.depwarn(
            """
            The `@mtkmodel` macro is deprecated. Please use the functional form with \
            `@components` instead.
            """, :mtkmodel
        )
        initial_conditions = __legacy_defaults__
    end

    if !(systems isa Vector{System})
        systems = Vector{System}(systems)
    end
    if !(eqs isa Vector{Equation})
        eqs = Equation[eqs]
    end
    eqs = eqs::Vector{Equation}

    iv = unwrap(iv)
    ps = vec(unwrap_vars(ps))
    dvs = vec(unwrap_vars(dvs))
    if iv !== nothing
        filter!(!Base.Fix2(_is_unknown_delay_or_evalat, iv), dvs)
    end
    brownians = unwrap_vars(brownians)

    if noise_eqs !== nothing
        noise_eqs = unwrap_vars(noise_eqs)
    end

    costs = vec(unwrap_vars(costs))

    if !(inputs isa OrderedSet{SymbolicT})
        inputs = unwrap.(inputs)
        inputs = OrderedSet{SymbolicT}(inputs)
    end
    if !(outputs isa OrderedSet{SymbolicT})
        outputs = unwrap.(outputs)
        outputs = OrderedSet{SymbolicT}(outputs)
    end
    for subsys in systems
        for var in get_inputs(subsys)
            push!(inputs, renamespace(subsys, var))
        end
        for var in get_outputs(subsys)
            push!(outputs, renamespace(subsys, var))
        end
    end
    var_to_name = Dict{Symbol, SymbolicT}()

    bindings = defsdict(bindings)
    initial_conditions = defsdict(initial_conditions)
    guesses = defsdict(guesses)
    all_dvs = as_atomic_array_set(dvs)
    if iv === nothing
        for k in keys(bindings)
            k in all_dvs || continue
            throw(
                ArgumentError(
                    """
                    Bindings for variables are enforced during initialization. Since \
                    time-independent systems only perform parameter initialization, \
                    bindings for variables in such systems are invalid. $k was found to have \
                    a binding in the system $name.
                    """
                )
            )
        end
    end
    let initial_conditions = discover_from_metadata ? initial_conditions : SymmapT(),
            bindings = discover_from_metadata ? bindings : SymmapT(),
            guesses = discover_from_metadata ? guesses : SymmapT(),
            inputs = discover_from_metadata ? inputs : OrderedSet{SymbolicT}(),
            outputs = discover_from_metadata ? outputs : OrderedSet{SymbolicT}()

        process_variables!(var_to_name, initial_conditions, bindings, guesses, dvs)
        process_variables!(var_to_name, initial_conditions, bindings, guesses, ps)
        buffer = SymbolicT[]
        for eq in observed
            push!(buffer, eq.lhs)
            if !(iv isa SymbolicT && _is_unknown_delay_or_evalat(eq.rhs, iv)) &&
                    !iscalledparameter(eq.rhs)
                push!(buffer, eq.rhs)
            end
        end
        process_variables!(var_to_name, initial_conditions, bindings, guesses, buffer)

        for var in dvs
            if isinput(var)
                push!(inputs, var)
            elseif isoutput(var)
                push!(outputs, var)
            end
        end
    end

    state_priorities = as_atomicmap(Int, state_priorities)
    irreducibles = parse_atomicset(irreducibles)
    if discover_from_metadata
        collect_metadata!(VariableStatePriority, state_priorities, dvs)
        collect_metadata!(VariableIrreducible, irreducibles, dvs)
    end

    filter!(!(Base.Fix1(===, COMMON_NOTHING) ∘ last), initial_conditions)
    filter!(!(Base.Fix1(===, COMMON_NOTHING) ∘ last), bindings)
    filter!(!(Base.Fix1(===, COMMON_NOTHING) ∘ last), guesses)

    if iv === nothing
        filterer = let initial_conditions = initial_conditions, all_dvs = all_dvs
            function _filterer(kvp)
                k = kvp[1]
                if k in all_dvs
                    initial_conditions[k] = kvp[2]
                    return false
                end
                return true
            end
        end
        filter!(filterer, bindings)
    end

    check_bindings(ps, bindings)
    bindings = ROSymmapT(bindings)

    if !allunique(map(nameof, systems))
        nonunique_subsystems(systems)
    end
    continuous_events,
        discrete_events = create_symbolic_events(
        continuous_events, discrete_events
    )

    if iv === nothing && (!isempty(continuous_events) || !isempty(discrete_events))
        throw(EventsInTimeIndependentSystemError(continuous_events, discrete_events))
    end

    if is_dde === nothing
        is_dde = _check_if_dde(eqs, iv, systems)
    end

    _assertions = Dict{SymbolicT, String}()
    for (k, v) in assertions
        _assertions[unwrap(k)::SymbolicT] = v
    end
    assertions = _assertions

    if isempty(metadata)
        metadata = MetadataT()
    elseif metadata isa MetadataT
        metadata = metadata
    else
        meta = MetadataT()
        for kvp in metadata
            meta = Base.ImmutableDict(meta, kvp)
        end
        metadata = meta
    end
    metadata = refreshed_metadata(metadata)
    jumps = Vector{JumpType}(jumps)
    return System(
        __get_new_tag(), eqs, noise_eqs, jumps, constraints,
        costs, consolidate, dvs, ps, brownians, iv, observed,
        var_to_name, name, description, bindings, initial_conditions, guesses, systems, initialization_eqs,
        continuous_events, discrete_events, connector_type, assertions, metadata, gui_metadata, is_dde,
        tstops, inputs, outputs, tearing_state, true, false,
        nothing, nothing, ignored_connections, preface, parent,
        initializesystem, is_initializesystem, is_discrete, state_priorities, irreducibles; checks
    )
end

function _is_unknown_delay_or_evalat(x::SymbolicT, iv::SymbolicT)
    x = split_indexed_var(x)[1]
    return Moshi.Match.@match x begin
        BSImpl.Term(; f, args) && if f isa SymbolicT end => begin
            !isequal(args[1], iv) && Moshi.Match.@match args[1] begin
                BSImpl.Term(; f = f2) => !(f2 isa SymbolicT) || SU.is_function_symbolic(f2)
                _ => true
            end
        end
        BSImpl.Term(; f, args) && if f === real || f === imag end => begin
            _is_unknown_delay_or_evalat(args[1], iv)
        end
        BSImpl.Term(; f, args) && if f isa Differential end => begin
            _is_unknown_delay_or_evalat(args[1], iv)
        end
        _ => false
    end
end

@noinline function nonunique_subsystems(systems)
    sysnames = nameof.(systems)
    unique_sysnames = Set(sysnames)
    throw(NonUniqueSubsystemsError(sysnames, unique_sysnames))
end

SymbolicIndexingInterface.getname(x::AbstractSystem) = nameof(x)

"""
    $(TYPEDSIGNATURES)

Create a time-independent [`System`](@ref) with the given equations `eqs`, unknowns `dvs`
and parameters `ps`.
"""
function System(eqs::Vector{Equation}, dvs, ps; kwargs...)
    return System(eqs, nothing, dvs, ps; kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Create a time-dependent system with the given equations `eqs` and independent variable `iv`.
Discover variables, parameters and brownians in the system by parsing the equations and
other symbolic expressions passed to the system.
"""
function System(eqs::Vector{Equation}, iv; kwargs...)
    iv === nothing && return System(eqs; kwargs...)

    diffvars = OrderedSet{SymbolicT}()
    othervars = OrderedSet{SymbolicT}()
    ps = OrderedSet{SymbolicT}()
    diffeqs = Equation[]
    othereqs = Equation[]
    iv = unwrap(iv)
    for eq in eqs
        collect_vars!(othervars, ps, eq, iv)
        if iscall(eq.lhs) && operation(eq.lhs) isa Differential
            var, _ = var_from_nested_derivative(eq.lhs)
            if var in diffvars
                throw(
                    ArgumentError(
                        """
                            The differential variable $var is not unique in the system of \
                        equations.
                        """
                    )
                )
            end
            # this check ensures var is correctly scoped, since `collect_vars!` won't pick
            # it up if it belongs to an ancestor system.
            if var in othervars
                push!(diffvars, var)
            end
            push!(diffeqs, eq)
        else
            push!(othereqs, eq)
        end
    end

    allunknowns = union(diffvars, othervars)
    eqs = [diffeqs; othereqs]

    brownians = Set{SymbolicT}()
    for x in allunknowns
        x = unwrap(x)
        if getvariabletype(x) == BROWNIAN
            push!(brownians, x)
        end
    end
    setdiff!(allunknowns, brownians)

    cstrs = Vector{Union{Equation, Inequality}}(get(kwargs, :constraints, []))
    _cstrunknowns, cstrps = process_constraint_system(cstrs, allunknowns, ps, iv)
    cstrunknowns = empty(_cstrunknowns)
    for var in _cstrunknowns
        op, inner = Moshi.Match.@match var begin
            BSImpl.Term(; f, args) && if f isa Differential end => (f, args[1])
            _ => begin
                push!(cstrunknowns, var)
                continue
            end
        end
        Moshi.Match.@match inner begin
            BSImpl.Term(; f, args) && if f isa SymbolicT && SU.isconst(args[1]) end => begin
                push!(cstrunknowns, op(f(iv)))
            end
            _ => push!(cstrunknowns, var)
        end
    end
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

    noiseeqs = get(kwargs, :noise_eqs, nothing)
    if noiseeqs !== nothing
        # validate noise equations
        noisedvs = OrderedSet{SymbolicT}()
        noiseps = OrderedSet{SymbolicT}()
        collect_vars!(noisedvs, noiseps, noiseeqs, iv)
        for dv in noisedvs
            dv ∈ allunknowns ||
                throw(ArgumentError("Variable $dv in noise equations is not an unknown of the system."))
        end
    end

    return System(
        eqs, iv, collect(allunknowns), collect(new_ps), collect(brownians); kwargs...
    )
end

"""
    $(TYPEDSIGNATURES)

Create a time-independent system with the given equations `eqs`. Discover variables and
parameters in the system by parsing the equations and other symbolic expressions passed to
the system.
"""
function System(eqs::Vector{Equation}; kwargs...)
    eqs = collect(eqs)

    allunknowns = OrderedSet{SymbolicT}()
    ps = OrderedSet{SymbolicT}()
    for eq in eqs
        collect_vars!(allunknowns, ps, eq, nothing)
    end
    for ssys in get(kwargs, :systems, System[])
        collect_scoped_vars!(allunknowns, ps, ssys, nothing)
    end
    costs = get(kwargs, :costs, [])
    for val in costs
        collect_vars!(allunknowns, ps, val, nothing)
    end

    cstrs = Vector{Union{Equation, Inequality}}(get(kwargs, :constraints, []))
    for eq in cstrs
        collect_vars!(allunknowns, ps, eq, nothing)
    end

    new_ps = gather_array_params(ps)

    return System(eqs, nothing, collect(allunknowns), collect(new_ps); kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Create a `System` with a single equation `eq`.
"""
System(eq::Equation, args...; kwargs...) = System([eq], args...; kwargs...)

function gather_array_params(ps)
    new_ps = OrderedSet{SymbolicT}()
    for p in ps
        arr, isarr = split_indexed_var(p)
        sh = SU.shape(arr)
        if isarr
            if !(sh isa SU.Unknown) && all(in(ps) ∘ Base.Fix1(getindex, arr), SU.stable_eachindex(arr))
                push!(new_ps, arr)
            else
                push!(new_ps, p)
            end
        else
            if sh isa SU.ShapeVecT && !isempty(sh)
                for i in SU.stable_eachindex(arr)
                    delete!(new_ps, arr[i])
                end
            end
            push!(new_ps, p)
        end
    end
    return new_ps
end

struct ConstraintValidator
    innervars::Set{SymbolicT}
end
ConstraintValidator() = ConstraintValidator(Set{SymbolicT}())

function (cv::ConstraintValidator)(x::SymbolicT)
    Moshi.Match.@match x begin
        BSImpl.Term(; f, args) && if f isa SymbolicT && !SU.is_function_symbolic(f) end => begin
            empty!(cv.innervars)
            SU.search_variables!(cv.innervars, args[1])
            for var in cv.innervars
                Moshi.Match.@match var begin
                    BSImpl.Term() => throw(
                        ArgumentError(
                            """
                            Arguments of a delayed or evaluated (via `EvalAt`) variable \
                            cannot be other dependent variables. Found $x which \
                            contains $var.
                            """
                        )
                    )
                    _ => nothing
                end
            end
        end
        _ => nothing
    end
    return false
end

"""
Process variables in constraints of the (ODE) System.
"""
function process_constraint_system(
        constraints::Vector{Union{Equation, Inequality}}, sts, ps, iv; validate = true,
    )
    isempty(constraints) && return OrderedSet{SymbolicT}(), OrderedSet{SymbolicT}()

    constraintsts = OrderedSet{SymbolicT}()
    constraintps = OrderedSet{SymbolicT}()
    validator = ConstraintValidator()
    for cons in constraints
        collect_vars!(constraintsts, constraintps, cons, iv)
        union!(constraintsts, collect_applied_operators(cons, Differential))
        SU.query(validator, cons.lhs)
        SU.query(validator, cons.rhs)
    end

    # Validate the states.
    if validate
        validate_vars_and_find_ps!(constraintsts, constraintps, sts, iv)
    end

    return constraintsts, constraintps
end

"""
Process the costs for the constraint system.
"""
function process_costs(costs::Vector, sts, ps, iv)
    coststs = OrderedSet{SymbolicT}()
    costps = OrderedSet{SymbolicT}()
    for cost in costs
        collect_vars!(coststs, costps, cost, iv)
    end

    validate_vars_and_find_ps!(coststs, costps, sts, iv)
    return coststs, costps
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
            SU.query(isequal(iv), var) && (
                var ∈ sts ||
                    throw(ArgumentError("Time-dependent variable $var is not an unknown of the system."))
            )
        elseif length(arguments(var)) > 1
            throw(ArgumentError("Too many arguments for variable $var."))
        elseif length(arguments(var)) == 1
            if iscall(var) && operation(var) isa Differential
                var = only(arguments(var))
            end
            arg = only(arguments(var))
            operation(var)(iv) ∈ sts ||
                throw(ArgumentError("Variable $var is not a variable of the System. Called variables must be variables of the System."))

            isequal(arg, iv) || isparameter(arg) || isconst(arg) && symtype(arg) <: Real ||
                throw(ArgumentError("Invalid argument specified for variable $var. The argument of the variable should be either $iv, a parameter, or a value specifying the time that the constraint holds."))

            isparameter(arg) && !isequal(arg, iv) && push!(auxps, arg)
        else
            var ∈ sts &&
                @warn "Variable $var has no argument. It will be interpreted as $var($iv), and the constraint will apply to the entire interval."
        end
    end
    return
end

"""
    $(TYPEDSIGNATURES)

Check if a system is a (possibly implicit) discrete system. Hybrid systems are turned into
callbacks, so checking if any LHS is shifted is sufficient. If a variable is shifted in
the input equations there _will_ be a `Shift` equation in the simplified system.
"""
function is_discrete_system(sys::System)
    return get_is_discrete(sys) || any(eq -> isoperator(eq.lhs, Shift), equations(sys))
end

SymbolicIndexingInterface.is_time_dependent(sys::System) = get_iv(sys) !== nothing

"""
    is_dde(sys::AbstractSystem)

Return a boolean indicating whether a system represents a set of delay
differential equations.
"""
is_dde(sys::AbstractSystem) = has_is_dde(sys) && get_is_dde(sys)

_check_if_dde(eqs::Vector{Equation}, iv::Nothing, subsystems::Vector{System}) = false
function _check_if_dde(eqs::Vector{Equation}, iv::SymbolicT, subsystems::Vector{System})
    any(ModelingToolkitBase.is_dde, subsystems) && return true
    pred = Base.Fix2(isdelay, iv)
    for eq in eqs
        SU.query(pred, eq.lhs) && return true
        SU.query(pred, eq.rhs) && return true
    end
    return false
end

"""
    $(TYPEDSIGNATURES)

Flatten the hierarchical structure of a system, collecting all equations, unknowns, etc.
into one top-level system after namespacing appropriately.
"""
function flatten(sys::System, noeqs = false)
    systems = get_systems(sys)
    isempty(systems) && return sys
    costs = cost(sys)
    if _iszero(costs)
        costs = SymbolicT[]
    else
        costs = SymbolicT[costs]
    end
    # We don't include `ignored_connections` in the flattened system, because
    # connection expansion inherently requires the hierarchy structure. If the system
    # is being flattened, then we no longer want to expand connections (or have already
    # done so) and thus don't care about `ignored_connections`.
    return System(
        noeqs ? Equation[] : equations(sys), get_iv(sys), unknowns(sys),
        parameters(sys; initial_parameters = true), brownians(sys);
        jumps = jumps(sys), constraints = constraints(sys), costs = costs,
        consolidate = default_consolidate, observed = observed(sys),
        bindings = bindings(sys), initial_conditions = initial_conditions(sys),
        guesses = guesses(sys),
        continuous_events = continuous_events(sys),
        discrete_events = discrete_events(sys), assertions = assertions(sys),
        is_dde = is_dde(sys), tstops = symbolic_tstops(sys),
        initialization_eqs = initialization_equations(sys),
        inputs = inputs(sys), outputs = outputs(sys),
        state_priorities = state_priorities(sys),
        irreducibles = irreducibles(sys),
        # without this, any initial conditions/bindings/guesses obtained from metadata that
        # were later removed by the user will be re-added. Right now, we just want to
        # retain `initial_conditions(sys)` as-is.
        discover_from_metadata = false, metadata = get_metadata(sys),
        gui_metadata = get_gui_metadata(sys),
        description = description(sys), name = nameof(sys)
    )
end

has_massactionjumps(js::System) = any(x -> x isa MassActionJump, jumps(js))
has_constantratejumps(js::System) = any(x -> x isa ConstantRateJump, jumps(js))
has_variableratejumps(js::System) = any(x -> x isa VariableRateJump, jumps(js))
# TODO: do we need this? it's kind of weird to keep
has_equations(js::System) = !isempty(equations(js))

function noise_equations_equal(sys1::System, sys2::System)
    neqs1 = get_noise_eqs(sys1)
    neqs2 = get_noise_eqs(sys2)
    if neqs1 === nothing && neqs2 === nothing
        return true
    elseif neqs1 === nothing || neqs2 === nothing
        return false
    end
    ndims(neqs1) == ndims(neqs2) || return false

    eqs1 = get_eqs(sys1)
    eqs2 = get_eqs(sys2)

    # get the permutation vector of `eqs2` in terms of `eqs1`
    # eqs1_used tracks the elements of `eqs1` already used in the permutation
    eqs1_used = falses(length(eqs1))
    # the permutation of `eqs1` that gives `eqs2`
    eqs2_perm = Int[]
    for eq in eqs2
        # find the first unused element of `eqs1` equal to `eq`
        idx = findfirst(i -> isequal(eq, eqs1[i]) && !eqs1_used[i], eachindex(eqs1))
        # none found, so noise equations are not equal
        idx === nothing && return false
        push!(eqs2_perm, idx)
    end

    if neqs1 isa Vector
        return isequal(@view(neqs1[eqs2_perm]), neqs2)
    else
        return isequal(@view(neqs1[eqs2_perm, :]), neqs2)
    end
end

function ignored_connections_equal(sys1::System, sys2::System)
    ic1 = get_ignored_connections(sys1)
    ic2 = get_ignored_connections(sys2)
    if ic1 === nothing && ic2 === nothing
        return true
    elseif ic1 === nothing || ic2 === nothing
        return false
    end
    return _eq_unordered(ic1[1], ic2[1]) && _eq_unordered(ic1[2], ic2[2])
end

"""
    $(TYPEDSIGNATURES)

Get the metadata associated with key `k` in system `sys` or `default` if it does not exist.
"""
function SymbolicUtils.getmetadata(sys::AbstractSystem, k::DataType, default)
    meta = get_metadata(sys)
    return get(meta, k, default)
end

"""
    $(TYPEDSIGNATURES)

Set the metadata associated with key `k` in system `sys` to value `v`. This is an
out-of-place operation, and will return a shallow copy of `sys` with the appropriate
metadata values.
"""
function SymbolicUtils.setmetadata(sys::AbstractSystem, k::DataType, v)
    meta = get_metadata(sys)
    meta = Base.ImmutableDict(meta, k => v)::MetadataT
    return @set sys.metadata = meta
end

function SymbolicUtils.hasmetadata(sys::AbstractSystem, k::DataType)
    meta = get_metadata(sys)
    return haskey(meta, k)
end

"""
    $(TYPEDSIGNATURES)

Metadata key for systems containing the `problem_type` to be passed to the problem
constructor, where applicable. For example, if `getmetadata(sys, ProblemTypeCtx, nothing)`
is `CustomType()` then `ODEProblem(sys, ...).problem_type` will be `CustomType()` instead
of `StandardODEProblem`.
"""
struct ProblemTypeCtx end

"""
    $(TYPEDSIGNATURES)
"""
function check_complete(sys::System, obj)
    return iscomplete(sys) || throw(SystemNotCompleteError(obj))
end

"""
    $(TYPEDSIGNATURES)

Given a time-dependent system `sys` of ODEs, convert it to a time-independent system of
nonlinear equations that solve for the steady-state of the unknowns. This is done by
replacing every derivative `D(x)` of an unknown `x` with zero. Note that this process
does not retain noise equations, brownian terms, jumps or costs associated with `sys`.
All other information such as initial conditions, bindings, guesses, observed and
initialization equations are retained. The independent variable of `sys` becomes a
parameter of the returned system.

If `sys` is hierarchical (it contains subsystems) this transformation will be applied
recursively to all subsystems. The output system will be marked as `complete` if and only
if the input system is also `complete`. This also retains the `split` flag passed to
`complete`.

See also: [`complete`](@ref).
"""
function NonlinearSystem(sys::System)
    if !is_time_dependent(sys)
        throw(ArgumentError("`NonlinearSystem` constructor expects a time-dependent `System`"))
    end
    eqs = equations(sys)
    obs = observed(sys)
    D = Differential(get_iv(sys))
    subrules = Dict([D(x) => 0.0 for x in unknowns(sys)])
    for var in brownians(sys)
        subrules[var] = 0.0
    end
    eqs = map(eqs) do eq
        substitute(eq, subrules)
    end
    new_ps = [parameters(sys); get_iv(sys)]
    if iscomplete(sys)
        append!(new_ps, collect(bound_parameters(sys)))
    end
    nsys = System(
        eqs, unknowns(sys), new_ps;
        bindings = merge(bindings(sys), Dict(get_iv(sys) => Inf)),
        initial_conditions = initial_conditions(sys), guesses = guesses(sys),
        initialization_eqs = initialization_equations(sys), name = nameof(sys),
        observed = obs, systems = map(NonlinearSystem, get_systems(sys))
    )
    if iscomplete(sys)
        nsys = complete(nsys; split = is_split(sys))
    end
    return nsys
end

########
# Utility constructors
########

"""
    $(TYPEDSIGNATURES)

Construct a time-independent [`System`](@ref) for optimizing the specified scalar `cost`.
The system will have no equations.

Unknowns and parameters of the system are inferred from the cost and other values (such as
initial conditions) passed to it.

All keyword arguments are the same as those of the [`System`](@ref) constructor.
"""
function OptimizationSystem(cost; kwargs...)
    return System(Equation[]; costs = [cost], kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Identical to the corresponding single-argument `OptimizationSystem` constructor, except
the unknowns and parameters are specified by passing arrays of symbolic variables to `dvs`
and `ps` respectively.
"""
function OptimizationSystem(cost, dvs, ps; kwargs...)
    return System(Equation[], nothing, dvs, ps; costs = [cost], kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Construct a time-independent [`System`](@ref) for optimizing the specified multi-objective
`cost`. The cost will be reduced to a scalar using the `consolidate` function. This
defaults to summing the specified cost and that of all subsystems. The system will have no
equations.

Unknowns and parameters of the system are inferred from the cost and other values (such as
initial conditions) passed to it.

All keyword arguments are the same as those of the [`System`](@ref) constructor.
"""
function OptimizationSystem(cost::Array; kwargs...)
    return System(Equation[]; costs = vec(cost), kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Identical to the corresponding single-argument `OptimizationSystem` constructor, except
the unknowns and parameters are specified by passing arrays of symbolic variables to `dvs`
and `ps` respectively.
"""
function OptimizationSystem(cost::Array, dvs, ps; kwargs...)
    return System(Equation[], nothing, dvs, ps; costs = vec(cost), kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Construct a [`System`](@ref) to solve a system of jump equations. `jumps` is an array of
jumps, expressed using `JumpProcesses.MassActionJump`, `JumpProcesses.ConstantRateJump`
and `JumpProcesses.VariableRateJump`. It can also include standard equations to simulate
jump-diffusion processes. `iv` should be the independent variable of the system.

All keyword arguments are the same as those of the [`System`](@ref) constructor.
"""
function JumpSystem(jumps, iv; kwargs...)
    mask = isa.(jumps, Equation)
    eqs = Vector{Equation}(jumps[mask])
    jumps = jumps[.!mask]
    return System(eqs, iv; jumps, kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Identical to the 2-argument `JumpSystem` constructor, but uses the explicitly provided
`dvs` and `ps` for unknowns and parameters of the system.
"""
function JumpSystem(jumps, iv, dvs, ps; kwargs...)
    mask = isa.(jumps, Equation)
    eqs = Vector{Equation}(jumps[mask])
    jumps = jumps[.!mask]
    return System(eqs, iv, dvs, ps; jumps, kwargs...)
end

"""
    SDESystem(eqs::Vector{Equation}, noise, iv; is_scalar_noise = false, kwargs...)

Construct a system of equations with associated noise terms. Instead of specifying noise
using [`@brownians`](@ref) variables, it is specified using a noise matrix `noise`. `iv` is
the independent variable of the system.

In the general case, `noise` should be a `N x M` matrix where `N` is the number of
equations (`length(eqs)`) and `M` is the number of independent random variables.
`noise[i, j]` is the diffusion term for equation `i` and random variable `j`. If the noise
is diagonal (`N == M` and `noise[i, j] == 0` for all `i != j`) it can be specified as a
`Vector` of length `N` corresponding to the diagonal of the noise matrix. As a special
case, if all equations have the same noise then all rows of `noise` are identical. This
is known as "scalar noise". In this case, `noise` can be a `Vector` corresponding to the
repeated row and `is_scalar_noise` must be `true`.

Note that systems created in this manner cannot be used hierarchically. This should only
be used to construct flattened systems. To use such a system hierarchically, it must be
converted to use brownian variables using [`noise_to_brownians`](@ref). [`mtkcompile`](@ref)
will automatically perform this conversion.

All keyword arguments are the same as those of the [`System`](@ref) constructor.
"""
function SDESystem(
        eqs::Vector{Equation}, noise, iv; is_scalar_noise = false,
        kwargs...
    )
    if is_scalar_noise
        if !(noise isa Vector)
            throw(ArgumentError("Expected noise to be a vector if `is_scalar_noise`"))
        end
        noise = repeat(reshape(noise, (1, :)), length(eqs))
    end
    return System(eqs, iv; noise_eqs = noise, kwargs...)
end

"""
    SDESystem(eqs::Vector{Equation}, noise, iv, dvs, ps; is_scalar_noise = false, kwargs...)


Identical to the 3-argument `SDESystem` constructor, but uses the explicitly provided
`dvs` and `ps` for unknowns and parameters of the system.
"""
function SDESystem(
        eqs::Vector{Equation}, noise, iv, dvs, ps; is_scalar_noise = false,
        kwargs...
    )
    if is_scalar_noise
        if !(noise isa Vector)
            throw(ArgumentError("Expected noise to be a vector if `is_scalar_noise`"))
        end
        noise = repeat(reshape(noise, (1, :)), length(eqs))
    end
    return System(eqs, iv, dvs, ps; noise_eqs = noise, kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Attach the given noise matrix `noise` to the system `sys`.
"""
function SDESystem(sys::System, noise; kwargs...)
    return SDESystem(equations(sys), noise, get_iv(sys); kwargs...)
end

struct SystemNotCompleteError <: Exception
    obj::Any
end

function Base.showerror(io::IO, err::SystemNotCompleteError)
    return print(
        io, """
        A completed system is required. Call `complete` or `mtkcompile` on the \
        system before creating a `$(err.obj)`.
        """
    )
end

struct IllFormedNoiseEquationsError <: Exception
    noise_eqs_rows::Int
    eqs_length::Int
end

function Base.showerror(io::IO, err::IllFormedNoiseEquationsError)
    return print(
        io, """
        Noise equations are ill-formed. The number of rows must match the number of drift \
        equations. `size(neqs, 1) == $(err.noise_eqs_rows) != length(eqs) == \
        $(err.eqs_length)`.
        """
    )
end

function NoNameError()
    return ArgumentError(
        """
        The `name` keyword must be provided. Please consider using the `@named` macro.
        """
    )
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
    return
end

struct EventsInTimeIndependentSystemError <: Exception
    cevents::Vector
    devents::Vector
end

function Base.showerror(io::IO, err::EventsInTimeIndependentSystemError)
    return println(
        io, """
        Events are not supported in time-independent systems. Provide an independent variable to \
        make the system time-dependent or remove the events.

        The following continuous events were provided:
        $(err.cevents)

        The following discrete events were provided:
        $(err.devents)
        """
    )
end

function supports_initialization(sys::System)
    return isempty(get_systems(sys)) && isempty(jumps(sys)) &&
        isempty(get_costs(sys)) && isempty(get_constraints(sys))
end

safe_eachrow(::Nothing) = nothing
safe_eachrow(x::AbstractArray) = eachrow(x)

safe_issetequal(::Nothing, ::Nothing) = true
safe_issetequal(::Nothing, x) = false
safe_issetequal(x, ::Nothing) = false
safe_issetequal(x, y) = issetequal(x, y)

"""
    $(TYPEDSIGNATURES)

Check if two systems are about equal, to the extent that ModelingToolkitBase.jl supports. Note
that if this returns `true`, the systems are not guaranteed to be exactly equivalent
(unless `sysa === sysb`) but are highly likely to represent a similar mathematical problem.
If this returns `false`, the systems are very likely to be different.
"""
function Base.isapprox(sysa::System, sysb::System)
    sysa === sysb && return true
    return nameof(sysa) == nameof(sysb) &&
        isequal(get_iv(sysa), get_iv(sysb)) &&
        issetequal(get_eqs(sysa), get_eqs(sysb)) &&
        safe_issetequal(
        safe_eachrow(get_noise_eqs(sysa)), safe_eachrow(get_noise_eqs(sysb))
    ) &&
        issetequal(get_jumps(sysa), get_jumps(sysb)) &&
        issetequal(get_constraints(sysa), get_constraints(sysb)) &&
        issetequal(get_costs(sysa), get_costs(sysb)) &&
        isequal(get_consolidate(sysa), get_consolidate(sysb)) &&
        issetequal(get_unknowns(sysa), get_unknowns(sysb)) &&
        issetequal(get_ps(sysa), get_ps(sysb)) &&
        issetequal(get_brownians(sysa), get_brownians(sysb)) &&
        issetequal(get_observed(sysa), get_observed(sysb)) &&
        isequal(get_description(sysa), get_description(sysb)) &&
        isequal(get_bindings(sysa), get_bindings(sysb)) &&
        isequal(get_initial_conditions(sysa), get_initial_conditions(sysb)) &&
        isequal(get_guesses(sysa), get_guesses(sysb)) &&
        issetequal(get_initialization_eqs(sysa), get_initialization_eqs(sysb)) &&
        issetequal(get_continuous_events(sysa), get_continuous_events(sysb)) &&
        issetequal(get_discrete_events(sysa), get_discrete_events(sysb)) &&
        isequal(get_connector_type(sysa), get_connector_type(sysb)) &&
        isequal(get_assertions(sysa), get_assertions(sysb)) &&
        isequal(get_metadata(sysa), get_metadata(sysb)) &&
        isequal(get_is_dde(sysa), get_is_dde(sysb)) &&
        issetequal(get_tstops(sysa), get_tstops(sysb)) &&
        issetequal(get_inputs(sysa), get_inputs(sysb)) &&
        issetequal(get_outputs(sysa), get_outputs(sysb)) &&
        safe_issetequal(get_ignored_connections(sysa), get_ignored_connections(sysb)) &&
        isequal(get_is_initializesystem(sysa), get_is_initializesystem(sysb)) &&
        isequal(get_is_discrete(sysa), get_is_discrete(sysb)) &&
        isequal(get_state_priorities(sysa), get_state_priorities(sysb)) &&
        isequal(get_isscheduled(sysa), get_isscheduled(sysb))
end

_maybe_copy(x) = applicable(copy, x) ? copy(x) : x

function Base.copy(sys::System)
    return System(
        __get_new_tag(), copy(get_eqs(sys)), _maybe_copy(get_noise_eqs(sys)), copy(get_jumps(sys)),
        copy(get_constraints(sys)), copy(get_costs(sys)), get_consolidate(sys),
        copy(get_unknowns(sys)), copy(get_ps(sys)), copy(get_brownians(sys)), get_iv(sys),
        copy(get_observed(sys)), copy(get_var_to_name(sys)), nameof(sys), get_description(sys),
        copy(get_bindings(sys)), copy(get_initial_conditions(sys)), copy(get_guesses(sys)),
        map(copy, get_systems(sys)), copy(get_initialization_eqs(sys)),
        copy(get_continuous_events(sys)), copy(get_discrete_events(sys)), get_connector_type(sys),
        copy(get_assertions(sys)), refreshed_metadata(get_metadata(sys)), get_gui_metadata(sys),
        get_is_dde(sys), copy(get_tstops(sys)), copy(get_inputs(sys)), copy(get_outputs(sys)),
        get_tearing_state(sys), does_namespacing(sys), false, get_index_cache(sys),
        get_parameter_bindings_graph(sys), _maybe_copy(get_ignored_connections(sys)),
        _maybe_copy(get_preface(sys)), _maybe_copy(get_parent(sys)),
        _maybe_copy(get_initializesystem(sys)), get_is_initializesystem(sys),
        get_is_discrete(sys), copy(get_state_priorities(sys)), copy(get_irreducibles(sys)),
        copy(get_isscheduled(sys)), _maybe_copy(get_schedule(sys)); checks = false
    )
end
