struct Schedule{V <: BipartiteGraphs.Matching}
    """
    Maximal matching of variables to equations calculated during structural simplification.
    """
    var_eq_matching::V
    """
    Mapping of `Differential`s of variables to corresponding derivative expressions.
    """
    dummy_sub::Dict{Any, Any}
end

const MetadataT = Base.ImmutableDict{DataType, Any}

struct System <: AbstractSystem
    tag::UInt
    eqs::Vector{Equation}
    # nothing - no noise
    # vector - diagonal noise
    # matrix - generic form
    # column matrix - scalar noise
    noise_eqs::Union{Nothing, AbstractVector, AbstractMatrix}
    jumps::Vector{JumpType}
    constraints::Vector{Union{Equation, Inequality}}
    costs::Vector{<:Union{BasicSymbolic, Real}}
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
    metadata::MetadataT
    gui_metadata::Any # ?
    is_dde::Bool
    tstops::Vector{Any}
    tearing_state::Any
    namespacing::Bool
    complete::Bool
    index_cache::Union{Nothing, IndexCache}
    ignored_connections::Union{
        Nothing, Tuple{Vector{IgnoredAnalysisPoint}, Vector{IgnoredAnalysisPoint}}}
    preface::Any
    parent::Union{Nothing, System}
    initializesystem::Union{Nothing, System}
    is_initializesystem::Bool
    isscheduled::Bool
    schedule::Union{Schedule, Nothing}

    function System(
            tag, eqs, noise_eqs, jumps, constraints, costs, consolidate, unknowns, ps,
            brownians, iv, observed, parameter_dependencies, var_to_name, name, description,
            defaults, guesses, systems, initialization_eqs, continuous_events, discrete_events,
            connector_type, assertions = Dict{BasicSymbolic, String}(),
            metadata = MetadataT(), gui_metadata = nothing,
            is_dde = false, tstops = [], tearing_state = nothing, namespacing = true,
            complete = false, index_cache = nothing, ignored_connections = nothing,
            preface = nothing, parent = nothing, initializesystem = nothing,
            is_initializesystem = false, isscheduled = false, schedule = nothing;
            checks::Union{Bool, Int} = true)
        if is_initializesystem && iv !== nothing
            throw(ArgumentError("""
            Expected initialization system to be time-independent. Found independent
            variable $iv.
            """))
        end
        jumps = Vector{JumpType}(jumps)
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
        new(tag, eqs, noise_eqs, jumps, constraints, costs,
            consolidate, unknowns, ps, brownians, iv,
            observed, parameter_dependencies, var_to_name, name, description, defaults,
            guesses, systems, initialization_eqs, continuous_events, discrete_events,
            connector_type, assertions, metadata, gui_metadata, is_dde,
            tstops, tearing_state, namespacing, complete, index_cache, ignored_connections,
            preface, parent, initializesystem, is_initializesystem, isscheduled, schedule)
    end
end

function default_consolidate(costs, subcosts)
    # `reduce` instead of `sum` because the rrule for `sum` doesn't
    # handle the `init` kwarg.
    return reduce(+, costs; init = 0.0) + reduce(+, subcosts; init = 0.0)
end

function System(eqs::Vector{Equation}, iv, dvs, ps, brownians = [];
        constraints = Union{Equation, Inequality}[], noise_eqs = nothing, jumps = [],
        costs = BasicSymbolic[], consolidate = default_consolidate,
        observed = Equation[], parameter_dependencies = Equation[], defaults = Dict(),
        guesses = Dict(), systems = System[], initialization_eqs = Equation[],
        continuous_events = SymbolicContinuousCallback[], discrete_events = SymbolicDiscreteCallback[],
        connector_type = nothing, assertions = Dict{BasicSymbolic, String}(),
        metadata = MetadataT(), gui_metadata = nothing,
        is_dde = nothing, tstops = [], tearing_state = nothing,
        ignored_connections = nothing, parent = nothing,
        description = "", name = nothing, discover_from_metadata = true,
        initializesystem = nothing, is_initializesystem = false, preface = [],
        checks = true)
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

    costs = unwrap.(costs)
    if isempty(costs)
        costs = Union{BasicSymbolic, Real}[]
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
    continuous_events, discrete_events = create_symbolic_events(
        continuous_events, discrete_events, eqs, iv)

    if iv === nothing && (!isempty(continuous_events) || !isempty(discrete_events))
        throw(EventsInTimeIndependentSystemError(continuous_events, discrete_events))
    end

    if is_dde === nothing
        is_dde = _check_if_dde(eqs, iv, systems)
    end

    assertions = Dict{BasicSymbolic, String}(unwrap(k) => v for (k, v) in assertions)

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
    System(Threads.atomic_add!(SYSTEM_COUNT, UInt(1)), eqs, noise_eqs, jumps, constraints,
        costs, consolidate, dvs, ps, brownians, iv, observed, parameter_dependencies,
        var_to_name, name, description, defaults, guesses, systems, initialization_eqs,
        continuous_events, discrete_events, connector_type, assertions, metadata, gui_metadata, is_dde,
        tstops, tearing_state, true, false, nothing, ignored_connections, preface, parent,
        initializesystem, is_initializesystem; checks)
end

function System(eqs::Vector{Equation}, dvs, ps; kwargs...)
    System(eqs, nothing, dvs, ps; kwargs...)
end

function System(eqs::Vector{Equation}, iv; kwargs...)
    iv === nothing && return System(eqs; kwargs...)

    diffvars = OrderedSet()
    othervars = OrderedSet()
    ps = Set()
    diffeqs = Equation[]
    othereqs = Equation[]
    for eq in eqs
        if !(eq.lhs isa Union{Symbolic, Number, AbstractArray})
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
        # validate noise equations
        noisedvs = OrderedSet()
        noiseps = OrderedSet()
        collect_vars!(noisedvs, noiseps, noiseeqs, iv)
        for dv in noisedvs
            dv ∈ allunknowns ||
                throw(ArgumentError("Variable $dv in noise equations is not an unknown of the system."))
        end
    end

    return System(
        eqs, iv, collect(allunknowns), collect(new_ps), collect(brownians); kwargs...)
end

function System(eqs::Vector{Equation}; kwargs...)
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

System(eq::Equation, args...; kwargs...) = System([eq], args...; kwargs...)

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
        constraints::Vector{Union{Equation, Inequality}}, sts, ps, iv; validate = true)
    isempty(constraints) && return Set(), Set()

    constraintsts = OrderedSet()
    constraintps = OrderedSet()
    for cons in constraints
        collect_vars!(constraintsts, constraintps, cons, iv)
        union!(constraintsts, collect_applied_operators(cons, Differential))
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
                throw(ArgumentError("Variable $var is not a variable of the System. Called variables must be variables of the System."))

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
    is_dde(sys::AbstractSystem)

Return a boolean indicating whether a system represents a set of delay
differential equations.
"""
is_dde(sys::AbstractSystem) = has_is_dde(sys) && get_is_dde(sys)

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
    costs = cost(sys)
    if _iszero(costs)
        costs = Union{Real, BasicSymbolic}[]
    else
        costs = [costs]
    end
    # We don't include `ignored_connections` in the flattened system, because
    # connection expansion inherently requires the hierarchy structure. If the system
    # is being flattened, then we no longer want to expand connections (or have already
    # done so) and thus don't care about `ignored_connections`.
    return System(noeqs ? Equation[] : equations(sys), get_iv(sys), unknowns(sys),
        parameters(sys; initial_parameters = true), brownians(sys);
        jumps = jumps(sys), constraints = constraints(sys), costs = costs,
        consolidate = default_consolidate, observed = observed(sys),
        parameter_dependencies = parameter_dependencies(sys), defaults = defaults(sys),
        guesses = guesses(sys), continuous_events = continuous_events(sys),
        discrete_events = discrete_events(sys), assertions = assertions(sys),
        is_dde = is_dde(sys), tstops = symbolic_tstops(sys),
        initialization_eqs = initialization_equations(sys),
        # without this, any defaults/guesses obtained from metadata that were
        # later removed by the user will be re-added. Right now, we just want to
        # retain `defaults(sys)` as-is.
        discover_from_metadata = false, metadata = get_metadata(sys),
        description = description(sys), name = nameof(sys))
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

function Base.:(==)(sys1::System, sys2::System)
    sys1 === sys2 && return true
    iv1 = get_iv(sys1)
    iv2 = get_iv(sys2)
    isequal(iv1, iv2) &&
        isequal(nameof(sys1), nameof(sys2)) &&
        _eq_unordered(get_eqs(sys1), get_eqs(sys2)) &&
        noise_equations_equal(sys1, sys2) &&
        _eq_unordered(get_jumps(sys1), get_jumps(sys2)) &&
        _eq_unordered(get_constraints(sys1), get_constraints(sys2)) &&
        _eq_unordered(get_costs(sys1), get_costs(sys2)) &&
        isequal(get_consolidate(sys1), get_consolidate(sys2)) &&
        _eq_unordered(get_unknowns(sys1), get_unknowns(sys2)) &&
        _eq_unordered(get_ps(sys1), get_ps(sys2)) &&
        _eq_unordered(get_brownians(sys1), get_brownians(sys2)) &&
        _eq_unordered(get_observed(sys1), get_observed(sys2)) &&
        _eq_unordered(get_parameter_dependencies(sys1), get_parameter_dependencies(sys2)) &&
        isequal(get_description(sys1), get_description(sys2)) &&
        isequal(get_defaults(sys1), get_defaults(sys2)) &&
        isequal(get_guesses(sys1), get_guesses(sys2)) &&
        _eq_unordered(get_initialization_eqs(sys1), get_initialization_eqs(sys2)) &&
        _eq_unordered(get_continuous_events(sys1), get_continuous_events(sys2)) &&
        _eq_unordered(get_discrete_events(sys1), get_discrete_events(sys2)) &&
        isequal(get_connector_type(sys1), get_connector_type(sys2)) &&
        isequal(get_assertions(sys1), get_assertions(sys2)) &&
        isequal(get_metadata(sys1), get_metadata(sys2)) &&
        isequal(get_gui_metadata(sys1), get_gui_metadata(sys2)) &&
        get_is_dde(sys1) == get_is_dde(sys2) &&
        _eq_unordered(get_tstops(sys1), get_tstops(sys2)) &&
        # not comparing tearing states because checking if they're equal up to ordering
        # is difficult
        getfield(sys1, :namespacing) == getfield(sys2, :namespacing) &&
        getfield(sys1, :complete) == getfield(sys2, :complete) &&
        ignored_connections_equal(sys1, sys2) &&
        get_parent(sys1) == get_parent(sys2) &&
        get_isscheduled(sys1) == get_isscheduled(sys2) &&
        all(s1 == s2 for (s1, s2) in zip(get_systems(sys1), get_systems(sys2)))
end

function Base.hash(sys::System, h::UInt)
    h = hash(nameof(sys), h)
    h = hash(get_iv(sys), h)
    # be considerate of things compared using `_eq_unordered` in `==`
    eqs = get_eqs(sys)
    eq_sortperm = sortperm(eqs; by = string)
    h = hash(@view(eqs[eq_sortperm]), h)
    neqs = get_noise_eqs(sys)
    if neqs === nothing
        h = hash(nothing, h)
    elseif neqs isa Vector
        h = hash(@view(neqs[eq_sortperm]), h)
    else
        h = hash(@view(neqs[eq_sortperm, :]), h)
    end
    h = hash(Set(get_jumps(sys)), h)
    h = hash(Set(get_constraints(sys)), h)
    h = hash(Set(get_costs(sys)), h)
    h = hash(get_consolidate(sys), h)
    h = hash(Set(get_unknowns(sys)), h)
    h = hash(Set(get_ps(sys)), h)
    h = hash(Set(get_brownians(sys)), h)
    h = hash(Set(get_observed(sys)), h)
    h = hash(Set(get_parameter_dependencies(sys)), h)
    h = hash(get_description(sys), h)
    h = hash(get_defaults(sys), h)
    h = hash(get_guesses(sys), h)
    h = hash(Set(get_initialization_eqs(sys)), h)
    h = hash(Set(get_continuous_events(sys)), h)
    h = hash(Set(get_discrete_events(sys)), h)
    h = hash(get_connector_type(sys), h)
    h = hash(get_assertions(sys), h)
    h = hash(get_metadata(sys), h)
    h = hash(get_gui_metadata(sys), h)
    h = hash(get_is_dde(sys), h)
    h = hash(Set(get_tstops(sys)), h)
    h = hash(Set(getfield(sys, :namespacing)), h)
    h = hash(Set(getfield(sys, :complete)), h)
    ics = get_ignored_connections(sys)
    if ics === nothing
        h = hash(ics, h)
    else
        h = hash(Set(ics[1]), hash(Set(ics[2]), h), h)
    end
    h = hash(get_parent(sys), h)
    h = hash(get_isscheduled(sys), h)
    for s in get_systems(sys)
        h = hash(s, h)
    end
    return h
end

function SymbolicUtils.getmetadata(sys::AbstractSystem, k::DataType, default)
    meta = get_metadata(sys)
    return get(meta, k, default)
end

function SymbolicUtils.setmetadata(sys::AbstractSystem, k::DataType, v)
    meta = get_metadata(sys)
    meta = Base.ImmutableDict(meta, k => v)::MetadataT
    @set sys.metadata = meta
end

"""
    $(TYPEDSIGNATURES)
"""
function check_complete(sys::System, obj)
    iscomplete(sys) || throw(SystemNotCompleteError(obj))
end

function NonlinearSystem(sys::System)
    if !is_time_dependent(sys)
        throw(ArgumentError("`NonlinearSystem` constructor expects a time-dependent `System`"))
    end
    eqs = equations(sys)
    obs = observed(sys)
    subrules = Dict([D(x) => 0.0 for x in unknowns(sys)])
    eqs = map(eqs) do eq
        fast_substitute(eq, subrules)
    end
    nsys = System(eqs, unknowns(sys), [parameters(sys); get_iv(sys)];
        parameter_dependencies = parameter_dependencies(sys),
        defaults = merge(defaults(sys), Dict(get_iv(sys) => Inf)), guesses = guesses(sys),
        initialization_eqs = initialization_equations(sys), name = nameof(sys),
        observed = obs)
    if iscomplete(sys)
        nsys = complete(nsys; split = is_split(sys))
    end
    return nsys
end

########
# Utility constructors
########

function OptimizationSystem(cost; kwargs...)
    return System(Equation[]; costs = [cost], kwargs...)
end

function OptimizationSystem(cost, dvs, ps; kwargs...)
    return System(Equation[], nothing, dvs, ps; costs = [cost], kwargs...)
end

function OptimizationSystem(cost::Array; kwargs...)
    return System(Equation[]; costs = vec(cost), kwargs...)
end

function OptimizationSystem(cost::Array, dvs, ps; kwargs...)
    return System(Equation[], nothing, dvs, ps; costs = vec(cost), kwargs...)
end

function JumpSystem(jumps, iv; kwargs...)
    mask = isa.(jumps, Equation)
    eqs = Vector{Equation}(jumps[mask])
    jumps = jumps[.!mask]
    return System(eqs, iv; jumps, kwargs...)
end

function JumpSystem(jumps, iv, dvs, ps; kwargs...)
    mask = isa.(jumps, Equation)
    eqs = Vector{Equation}(jumps[mask])
    jumps = jumps[.!mask]
    return System(eqs, iv, dvs, ps; jumps, kwargs...)
end

function SDESystem(eqs::Vector{Equation}, noise, iv; is_scalar_noise = false, kwargs...)
    if is_scalar_noise
        if !(noise isa Vector)
            throw(ArgumentError("Expected noise to be a vector if `is_scalar_noise`"))
        end
        noise = repeat(reshape(noise, (1, :)), length(eqs))
    end
    return System(eqs, iv; noise_eqs = noise, kwargs...)
end

function SDESystem(
        eqs::Vector{Equation}, noise, iv, dvs, ps; is_scalar_noise = false, kwargs...)
    if is_scalar_noise
        if !(noise isa Vector)
            throw(ArgumentError("Expected noise to be a vector if `is_scalar_noise`"))
        end
        noise = repeat(reshape(noise, (1, :)), length(eqs))
    end
    return System(eqs, iv, dvs, ps; noise_eqs = noise, kwargs...)
end

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
