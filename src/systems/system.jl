struct Schedule
    var_sccs::Vector{Vector{Int}}
    dummy_sub::Dict{SymbolicT, SymbolicT}
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
        if is_initializesystem && iv !== nothing
            throw(ArgumentError())
        end
        @assert iv === nothing || symtype(iv) === Real
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
end
SymbolicIndexingInterface.getname(x::System) = nameof(x)
function System(eqs::Vector{Equation}, dvs, ps; kwargs...)
    System(eqs, nothing, dvs, ps; kwargs...)
end
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
                throw(ArgumentError())
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
function process_costs(costs::Vector, sts, ps, iv)
    coststs = OrderedSet{SymbolicT}()
    costps = OrderedSet{SymbolicT}()
    for cost in costs
        collect_vars!(coststs, costps, cost, iv)
    end
    validate_vars_and_find_ps!(coststs, costps, sts, iv)
    coststs, costps
end
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
function is_discrete_system(sys::System)
    get_is_discrete(sys) || any(eq -> isoperator(eq.lhs, Shift), equations(sys))
end
SymbolicIndexingInterface.is_time_dependent(sys::System) = get_iv(sys) !== nothing
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
function SymbolicUtils.getmetadata(sys::AbstractSystem, k::DataType, default)
    meta = get_metadata(sys)
    return get(meta, k, default)
end
function SymbolicUtils.setmetadata(sys::AbstractSystem, k::DataType, v)
    meta = get_metadata(sys)
    meta = Base.ImmutableDict(meta, k => v)::MetadataT
    @set sys.metadata = meta
end
function SymbolicUtils.hasmetadata(sys::AbstractSystem, k::DataType)
    meta = get_metadata(sys)
    haskey(meta, k)
end
struct ProblemTypeCtx end
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
function SDESystem(sys::System, noise; kwargs...)
    SDESystem(equations(sys), noise, get_iv(sys); kwargs...)
end
struct SystemNotCompleteError <: Exception
    obj::Any
end
function Base.showerror(io::IO, err::SystemNotCompleteError)
    print(io, )
end
struct IllFormedNoiseEquationsError <: Exception
    noise_eqs_rows::Int
    eqs_length::Int
end
function Base.showerror(io::IO, err::IllFormedNoiseEquationsError)
    print(io, )
end
function NoNameError()
    ArgumentError()
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
        io, )
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
