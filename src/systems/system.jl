struct Schedule
    var_sccs::Vector{Vector{Int}}
    """
    Mapping of `Differential`s of variables to corresponding derivative expressions.
    """
    dummy_sub::Dict{SymbolicT, SymbolicT}
end
const MetadataT = Base.ImmutableDict{DataType, Any}
abstract type MutableCacheKey end
const MutableCacheT = Dict{DataType, Any}
""""""
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
    All the explicit equations relating parameters. Equations here only contain parameters
    and are in the same format as `observed`.
    """
    parameter_dependencies::Vector{Equation}
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
    Default values that variables (unknowns/observables/parameters) should take when
    constructing a numerical problem from the system. These values can be overridden
    by initial values provided to the problem constructor. Defaults of parent systems
    take priority over those in child systems.
    """
    defaults::SymmapT
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
    gui_metadata::Any
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
        if is_initializesystem && iv !== nothing
            throw(ArgumentError("""
            Expected initialization system to be time-independent. Found independent
            variable $iv.
            """))
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
_sum_costs(costs::Vector{SymbolicT}) = SU.add_worker(VartypeT, costs)
_sum_costs(costs::Vector{Num}) = SU.add_worker(VartypeT, costs)
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
""""""
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
    name === nothing && throw(NoNameError())
    if !(systems isa Vector{System})
        systems = Vector{System}(systems)
    end
    if !(eqs isa Vector{Equation})
        eqs = Equation[eqs]
    end
    eqs = eqs::Vector{Equation}
    if !isempty(parameter_dependencies)
        @invokelatest warn_pdeps()
        append!(eqs, parameter_dependencies)
    end
    iv = unwrap(iv)
    ps = vec(unwrap_vars(ps))
    dvs = vec(unwrap_vars(dvs))
    if iv !== nothing
        filter!(!Base.Fix2(isdelay, iv), dvs)
    end
    brownians = unwrap_vars(brownians)
    if noise_eqs !== nothing
        noise_eqs = unwrap_vars(noise_eqs)
    end
    costs = vec(unwrap_vars(costs))
    defaults = defsdict(defaults)
    guesses = defsdict(guesses)
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
    let defaults = discover_from_metadata ? defaults : SymmapT(),
        guesses = discover_from_metadata ? guesses : SymmapT(),
        inputs = discover_from_metadata ? inputs : OrderedSet{SymbolicT}(),
        outputs = discover_from_metadata ? outputs : OrderedSet{SymbolicT}()
        process_variables!(var_to_name, defaults, guesses, dvs)
        process_variables!(var_to_name, defaults, guesses, ps)
        buffer = SymbolicT[]
        for eq in observed
            push!(buffer, eq.lhs)
            push!(buffer, eq.rhs)
        end
        process_variables!(var_to_name, defaults, guesses, buffer)
        for var in dvs
            if isinput(var)
                push!(inputs, var)
            elseif isoutput(var)
                push!(outputs, var)
            end
        end
    end
    filter!(!(Base.Fix1(===, COMMON_NOTHING) ∘ last), defaults)
    filter!(!(Base.Fix1(===, COMMON_NOTHING) ∘ last), guesses)
    if !allunique(map(nameof, systems))
        nonunique_subsystems(systems)
    end
    continuous_events,
    discrete_events = create_symbolic_events(
        continuous_events, discrete_events)
    if iv === nothing && (!isempty(continuous_events) || !isempty(discrete_events))
        throw(EventsInTimeIndependentSystemError(continuous_events, discrete_events))
    end
    if is_dde === nothing
        is_dde = _check_if_dde(eqs, iv, systems)
    end
    _assertions = Dict{SymbolicT, String}
    for (k, v) in assertions
        _assertions[unwrap(k)::SymbolicT] = v
    end
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
    System(Threads.atomic_add!(SYSTEM_COUNT, UInt(1)), eqs, noise_eqs, jumps, constraints,
        costs, consolidate, dvs, ps, brownians, iv, observed, Equation[],
        var_to_name, name, description, defaults, guesses, systems, initialization_eqs,
        continuous_events, discrete_events, connector_type, assertions, metadata, gui_metadata, is_dde,
        tstops, inputs, outputs, tearing_state, true, false,
        nothing, ignored_connections, preface, parent,
        initializesystem, is_initializesystem, is_discrete; checks)
end
@noinline function nonunique_subsystems(systems)
    sysnames = nameof.(systems)
    unique_sysnames = Set(sysnames)
    throw(NonUniqueSubsystemsError(sysnames, unique_sysnames))
end
@noinline function warn_pdeps()
    @warn """
    The `parameter_dependencies` keyword argument is deprecated. Please provide all
    such equations as part of the normal equations of the system.
    """
end
SymbolicIndexingInterface.getname(x::System) = nameof(x)
""""""
function System(eqs::Vector{Equation}, dvs, ps; kwargs...)
    System(eqs, nothing, dvs, ps; kwargs...)
end
""""""
function System(eqs::Vector{Equation}, iv; kwargs...)
    iv === nothing && return System(eqs; kwargs...)
    diffvars = OrderedSet{SymbolicT}()
    othervars = OrderedSet{SymbolicT}()
    ps = OrderedSet{SymbolicT}()
    diffeqs = Equation[]
    othereqs = Equation[]
    iv = unwrap(iv)
    for eq in eqs
        if !(eq.lhs isa Union{SymbolicT, Number, AbstractArray})
            push!(othereqs, eq)
            continue
        end
        collect_vars!(othervars, ps, eq, iv)
        if iscall(eq.lhs) && operation(eq.lhs) isa Differential
            var, _ = var_from_nested_derivative(eq.lhs)
            if var in diffvars
                throw(ArgumentError("""
                    The differential variable $var is not unique in the system of \
                    equations.
                """))
            end
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
    for eq in get(kwargs, :parameter_dependencies, Equation[])
        collect_vars!(allunknowns, ps, eq, iv)
    end
    cstrs = Vector{Union{Equation, Inequality}}(get(kwargs, :constraints, []))
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
    noiseeqs = get(kwargs, :noise_eqs, nothing)
    if noiseeqs !== nothing
        noisedvs = OrderedSet{SymbolicT}()
        noiseps = OrderedSet{SymbolicT}()
        collect_vars!(noisedvs, noiseps, noiseeqs, iv)
        for dv in noisedvs
            dv ∈ allunknowns ||
                throw(ArgumentError("Variable $dv in noise equations is not an unknown of the system."))
        end
    end
    return System(
        eqs, iv, collect(allunknowns), collect(new_ps), collect(brownians); kwargs...)
end
""""""
function System(eqs::Vector{Equation}; kwargs...)
    eqs = collect(eqs)
    allunknowns = OrderedSet{SymbolicT}()
    ps = OrderedSet{SymbolicT}()
    for eq in eqs
        collect_vars!(allunknowns, ps, eq, nothing)
    end
    for eq in get(kwargs, :parameter_dependencies, Equation[])
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
""""""
System(eq::Equation, args...; kwargs...) = System([eq], args...; kwargs...)
function gather_array_params(ps)
    new_ps = OrderedSet()
    for p in ps
        if iscall(p) && operation(p) === getindex
            par = arguments(p)[begin]
            if symbolic_has_known_size(p) && all(par[i] in ps for i in eachindex(par))
                push!(new_ps, par)
            else
                push!(new_ps, p)
            end
        else
            if symbolic_type(p) == ArraySymbolic() && symbolic_has_known_size(p)
                for i in eachindex(p)
                    delete!(new_ps, p[i])
                end
            end
            push!(new_ps, p)
        end
    end
    return new_ps
end
""""""
function process_constraint_system(
        constraints::Vector{Union{Equation, Inequality}}, sts, ps, iv; validate = true)
    isempty(constraints) && return OrderedSet{SymbolicT}(), OrderedSet{SymbolicT}()
    constraintsts = OrderedSet{SymbolicT}()
    constraintps = OrderedSet{SymbolicT}()
    for cons in constraints
        collect_vars!(constraintsts, constraintps, cons, iv)
        union!(constraintsts, collect_applied_operators(cons, Differential))
    end
    if validate
        validate_vars_and_find_ps!(constraintsts, constraintps, sts, iv)
    end
    return constraintsts, constraintps
end
""""""
function process_costs(costs::Vector, sts, ps, iv)
    coststs = OrderedSet{SymbolicT}()
    costps = OrderedSet{SymbolicT}()
    for cost in costs
        collect_vars!(coststs, costps, cost, iv)
    end
    validate_vars_and_find_ps!(coststs, costps, sts, iv)
    coststs, costps
end
""""""
function validate_vars_and_find_ps!(auxvars, auxps, sysvars, iv)
    sts = sysvars
    for var in auxvars
        if !iscall(var)
            SU.query(isequal(iv), var) && (var ∈ sts ||
             throw(ArgumentError("Time-dependent variable $var is not an unknown of the system.")))
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
end
""""""
function is_discrete_system(sys::System)
    get_is_discrete(sys) || any(eq -> isoperator(eq.lhs, Shift), equations(sys))
end
SymbolicIndexingInterface.is_time_dependent(sys::System) = get_iv(sys) !== nothing
""""""
is_dde(sys::AbstractSystem) = has_is_dde(sys) && get_is_dde(sys)
_check_if_dde(eqs::Vector{Equation}, iv::Nothing, subsystems::Vector{System}) = false
function _check_if_dde(eqs::Vector{Equation}, iv::SymbolicT, subsystems::Vector{System})
    any(ModelingToolkit.is_dde, subsystems) && return true
    pred = Base.Fix2(isdelay, iv)
    for eq in eqs
        SU.query(pred, eq.lhs) && return true
        SU.query(pred, eq.rhs) && return true
    end
    return false
end
""""""
function flatten(sys::System, noeqs = false)
    systems = get_systems(sys)
    isempty(systems) && return sys
    costs = cost(sys)
    if _iszero(costs)
        costs = SymbolicT[]
    else
        costs = SymbolicT[costs]
    end
    return System(noeqs ? Equation[] : equations(sys), get_iv(sys), unknowns(sys),
        parameters(sys; initial_parameters = true), brownians(sys);
        jumps = jumps(sys), constraints = constraints(sys), costs = costs,
        consolidate = default_consolidate, observed = observed(sys),
        defaults = defaults(sys), guesses = guesses(sys),
        continuous_events = continuous_events(sys),
        discrete_events = discrete_events(sys), assertions = assertions(sys),
        is_dde = is_dde(sys), tstops = symbolic_tstops(sys),
        initialization_eqs = initialization_equations(sys),
        inputs = inputs(sys), outputs = outputs(sys),
        discover_from_metadata = false, metadata = get_metadata(sys),
        description = description(sys), name = nameof(sys))
end
has_massactionjumps(js::System) = any(x -> x isa MassActionJump, jumps(js))
has_constantratejumps(js::System) = any(x -> x isa ConstantRateJump, jumps(js))
has_variableratejumps(js::System) = any(x -> x isa VariableRateJump, jumps(js))
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
    eqs1_used = falses(length(eqs1))
    eqs2_perm = Int[]
    for eq in eqs2
        idx = findfirst(i -> isequal(eq, eqs1[i]) && !eqs1_used[i], eachindex(eqs1))
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
""""""
function SymbolicUtils.getmetadata(sys::AbstractSystem, k::DataType, default)
    meta = get_metadata(sys)
    return get(meta, k, default)
end
""""""
function SymbolicUtils.setmetadata(sys::AbstractSystem, k::DataType, v)
    meta = get_metadata(sys)
    meta = Base.ImmutableDict(meta, k => v)::MetadataT
    @set sys.metadata = meta
end
function SymbolicUtils.hasmetadata(sys::AbstractSystem, k::DataType)
    meta = get_metadata(sys)
    haskey(meta, k)
end
""""""
struct ProblemTypeCtx end
""""""
function check_complete(sys::System, obj)
    iscomplete(sys) || throw(SystemNotCompleteError(obj))
end
""""""
function NonlinearSystem(sys::System)
    if !is_time_dependent(sys)
        throw(ArgumentError("`NonlinearSystem` constructor expects a time-dependent `System`"))
    end
    eqs = equations(sys)
    obs = observed(sys)
    subrules = Dict([D(x) => 0.0 for x in unknowns(sys)])
    for var in brownians(sys)
        subrules[var] = 0.0
    end
    eqs = map(eqs) do eq
        substitute(eq, subrules)
    end
    nsys = System(eqs, unknowns(sys), [parameters(sys); get_iv(sys)];
        defaults = merge(defaults(sys), Dict(get_iv(sys) => Inf)), guesses = guesses(sys),
        initialization_eqs = initialization_equations(sys), name = nameof(sys),
        observed = obs, systems = map(NonlinearSystem, get_systems(sys)))
    if iscomplete(sys)
        nsys = complete(nsys; split = is_split(sys))
        @set! nsys.parameter_dependencies = get_parameter_dependencies(sys)
    end
    return nsys
end
""""""
function OptimizationSystem(cost; kwargs...)
    return System(Equation[]; costs = [cost], kwargs...)
end
""""""
function OptimizationSystem(cost, dvs, ps; kwargs...)
    return System(Equation[], nothing, dvs, ps; costs = [cost], kwargs...)
end
""""""
function OptimizationSystem(cost::Array; kwargs...)
    return System(Equation[]; costs = vec(cost), kwargs...)
end
""""""
function OptimizationSystem(cost::Array, dvs, ps; kwargs...)
    return System(Equation[], nothing, dvs, ps; costs = vec(cost), kwargs...)
end
""""""
function JumpSystem(jumps, iv; kwargs...)
    mask = isa.(jumps, Equation)
    eqs = Vector{Equation}(jumps[mask])
    jumps = jumps[.!mask]
    return System(eqs, iv; jumps, kwargs...)
end
""""""
function JumpSystem(jumps, iv, dvs, ps; kwargs...)
    mask = isa.(jumps, Equation)
    eqs = Vector{Equation}(jumps[mask])
    jumps = jumps[.!mask]
    return System(eqs, iv, dvs, ps; jumps, kwargs...)
end
""""""
function SDESystem(eqs::Vector{Equation}, noise, iv; is_scalar_noise = false,
        parameter_dependencies = Equation[], kwargs...)
    if is_scalar_noise
        if !(noise isa Vector)
            throw(ArgumentError("Expected noise to be a vector if `is_scalar_noise`"))
        end
        noise = repeat(reshape(noise, (1, :)), length(eqs))
    end
    sys = System(eqs, iv; noise_eqs = noise, kwargs...)
    @set sys.parameter_dependencies = parameter_dependencies
end
""""""
function SDESystem(
        eqs::Vector{Equation}, noise, iv, dvs, ps; is_scalar_noise = false,
        parameter_dependencies = Equation[], kwargs...)
    if is_scalar_noise
        if !(noise isa Vector)
            throw(ArgumentError("Expected noise to be a vector if `is_scalar_noise`"))
        end
        noise = repeat(reshape(noise, (1, :)), length(eqs))
    end
    sys = System(eqs, iv, dvs, ps; noise_eqs = noise, kwargs...)
    @set sys.parameter_dependencies = parameter_dependencies
end
""""""
function SDESystem(sys::System, noise; kwargs...)
    SDESystem(equations(sys), noise, get_iv(sys); kwargs...)
end
struct SystemNotCompleteError <: Exception
    obj::Any
end
function Base.showerror(io::IO, err::SystemNotCompleteError)
    print(io, """
    A completed system is required. Call `complete` or `mtkcompile` on the \
    system before creating a `$(err.obj)`.
    """)
end
struct IllFormedNoiseEquationsError <: Exception
    noise_eqs_rows::Int
    eqs_length::Int
end
function Base.showerror(io::IO, err::IllFormedNoiseEquationsError)
    print(io, """
    Noise equations are ill-formed. The number of rows must match the number of drift \
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
    println(
        io, """
Events are not supported in time-independent systems. Provide an independent variable to \
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
safe_eachrow(::Nothing) = nothing
safe_eachrow(x::AbstractArray) = eachrow(x)
safe_issetequal(::Nothing, ::Nothing) = true
safe_issetequal(::Nothing, x) = false
safe_issetequal(x, ::Nothing) = false
safe_issetequal(x, y) = issetequal(x, y)
""""""
function Base.isapprox(sysa::System, sysb::System)
    sysa === sysb && return true
    return nameof(sysa) == nameof(sysb) &&
           isequal(get_iv(sysa), get_iv(sysb)) &&
           issetequal(get_eqs(sysa), get_eqs(sysb)) &&
           safe_issetequal(
               safe_eachrow(get_noise_eqs(sysa)), safe_eachrow(get_noise_eqs(sysb))) &&
           issetequal(get_jumps(sysa), get_jumps(sysb)) &&
           issetequal(get_constraints(sysa), get_constraints(sysb)) &&
           issetequal(get_costs(sysa), get_costs(sysb)) &&
           isequal(get_consolidate(sysa), get_consolidate(sysb)) &&
           issetequal(get_unknowns(sysa), get_unknowns(sysb)) &&
           issetequal(get_ps(sysa), get_ps(sysb)) &&
           issetequal(get_brownians(sysa), get_brownians(sysb)) &&
           issetequal(get_observed(sysa), get_observed(sysb)) &&
           issetequal(get_parameter_dependencies(sysa), get_parameter_dependencies(sysb)) &&
           isequal(get_description(sysa), get_description(sysb)) &&
           isequal(get_defaults(sysa), get_defaults(sysb)) &&
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
           isequal(get_isscheduled(sysa), get_isscheduled(sysb))
end
