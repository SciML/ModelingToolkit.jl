"""
$(TYPEDEF)

A system of ordinary differential equations.

# Fields
$(FIELDS)

# Example

```julia
using ModelingToolkit

@parameters σ ρ β
@variables t x(t) y(t) z(t)
D = Differential(t)

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]

@named de = ODESystem(eqs,t,[x,y,z],[σ,ρ,β],tspan=(0, 1000.0))
```
"""
struct ODESystem <: AbstractODESystem
    """
    tag: a tag for the system. If two systems have the same tag, then they are
    structurally ide, sym_to_stringntical.
    """
    tag::UInt
    """The ODEs defining the system."""
    eqs::Vector{Equation}
    """Independent variable."""
    iv::BasicSymbolic{Real}
    """
    Dependent (state) variables. Must not contain the independent variable.

    N.B.: If torn_matching !== nothing, this includes all variables. Actual
    ODE states are determined by the SelectedState() entries in `torn_matching`.
    """
    states::Vector
    """Parameter variables. Must not contain the independent variable."""
    ps::Vector
    """Time span."""
    tspan::Union{NTuple{2, Any}, Nothing}
    """Array variables."""
    var_to_name::Any
    """Control parameters (some subset of `ps`)."""
    ctrls::Vector
    """Observed states."""
    observed::Vector{Equation}
    """
    Time-derivative matrix. Note: this field will not be defined until
    [`calculate_tgrad`](@ref) is called on the system.
    """
    tgrad::RefValue{Vector{Num}}
    """
    Jacobian matrix. Note: this field will not be defined until
    [`calculate_jacobian`](@ref) is called on the system.
    """
    jac::RefValue{Any}
    """
    Control Jacobian matrix. Note: this field will not be defined until
    [`calculate_control_jacobian`](@ref) is called on the system.
    """
    ctrl_jac::RefValue{Any}
    """
    `Wfact` matrix. Note: this field will not be defined until
    [`generate_factorized_W`](@ref) is called on the system.
    """
    Wfact::RefValue{Matrix{Num}}
    """
    `Wfact_t` matrix. Note: this field will not be defined until
    [`generate_factorized_W`](@ref) is called on the system.
    """
    Wfact_t::RefValue{Matrix{Num}}
    """
    Name: the name of the system
    """
    name::Symbol
    """
    systems: The internal systems. These are required to have unique names.
    """
    systems::Vector{ODESystem}
    """
    defaults: The default values to use when initial conditions and/or
    parameters are not supplied in `ODEProblem`.
    """
    defaults::Dict
    """
    torn_matching: Tearing result specifying how to solve the system.
    """
    torn_matching::Union{Matching, Nothing}
    """
    connector_type: type of the system
    """
    connector_type::Any
    """
    preface: inject assignment statements before the evaluation of the RHS function.
    """
    preface::Any
    """
    continuous_events: A `Vector{SymbolicContinuousCallback}` that model events.
    The integrator will use root finding to guarantee that it steps at each zero crossing.
    """
    continuous_events::Vector{SymbolicContinuousCallback}
    """
    discrete_events: A `Vector{SymbolicDiscreteCallback}` that models events. Symbolic
    analog to `SciMLBase.DiscreteCallback` that executes an affect when a given condition is
    true at the end of an integration step.
    """
    discrete_events::Vector{SymbolicDiscreteCallback}
    """
    metadata: metadata for the system, to be used by downstream packages.
    """
    metadata::Any
    """
    gui_metadata: metadata for MTK GUI.
    """
    gui_metadata::Union{Nothing, GUIMetadata}
    """
    tearing_state: cache for intermediate tearing state
    """
    tearing_state::Any
    """
    substitutions: substitutions generated by tearing.
    """
    substitutions::Any
    """
    complete: if a model `sys` is complete, then `sys.x` no longer performs namespacing.
    """
    complete::Bool
    """
    discrete_subsystems: a list of discrete subsystems.
    """
    discrete_subsystems::Any
    """
    unknown_states: a list of actual states needed to be solved by solvers. Only
    used for ODAEProblem.
    """
    unknown_states::Union{Nothing, Vector{Any}}

    """
    maps states to indices in the parameter vector.
    """
    sym_to_index::Dict{Num, Int}

    function ODESystem(tag, deqs, iv, dvs, ps, tspan, var_to_name, ctrls, observed, tgrad,
        jac, ctrl_jac, Wfact, Wfact_t, name, systems, defaults,
        torn_matching, connector_type, preface, cevents,
        devents, metadata = nothing, gui_metadata = nothing,
        tearing_state = nothing,
        substitutions = nothing, complete = false,
        discrete_subsystems = nothing, unknown_states = nothing,
        sym_to_index = Dict{Num, Int}();
        checks::Union{Bool, Int} = true)
        if checks == true || (checks & CheckComponents) > 0
            check_variables(dvs, iv)
            check_parameters(ps, iv)
            check_equations(deqs, iv)
            check_equations(equations(cevents), iv)
        end
        if checks == true || (checks & CheckUnits) > 0
            all_dimensionless([dvs; ps; iv]) || check_units(deqs)
        end

        self = new(tag, deqs, iv, dvs, ps, tspan, var_to_name, ctrls, observed, tgrad, jac,
            ctrl_jac, Wfact, Wfact_t, name, systems, defaults, torn_matching,
            connector_type, preface, cevents, devents, metadata, gui_metadata,
            tearing_state, substitutions, complete, discrete_subsystems,
            unknown_states, sym_to_index)
        fill_unknown_states!(self)
    end
end

function fill_unknown_states!(sys::ODESystem)
    us = unknown_states(sys)
    sym_to_index = sys.sym_to_index
    sizehint!(sym_to_index, length(us))
    for (i, p) in enumerate(us)
        sym_to_index[p] = i
    end
    return sys
end

function ODESystem(deqs::AbstractVector{<:Equation}, iv, dvs, ps;
    controls = Num[],
    observed = Equation[],
    systems = ODESystem[],
    tspan = nothing,
    name = nothing,
    default_u0 = Dict(),
    default_p = Dict(),
    defaults = _merge(Dict(default_u0), Dict(default_p)),
    connector_type = nothing,
    preface = nothing,
    continuous_events = nothing,
    discrete_events = nothing,
    checks = true,
    metadata = nothing,
    gui_metadata = nothing)
    name === nothing &&
        throw(ArgumentError("The `name` keyword must be provided. Please consider using the `@named` macro"))
    deqs = scalarize(deqs)
    @assert all(control -> any(isequal.(control, ps)), controls) "All controls must also be parameters."

    iv′ = value(scalarize(iv))
    dvs′ = value.(scalarize(dvs))
    ps′ = value.(scalarize(ps))
    ctrl′ = value.(scalarize(controls))

    if !(isempty(default_u0) && isempty(default_p))
        Base.depwarn("`default_u0` and `default_p` are deprecated. Use `defaults` instead.",
            :ODESystem, force = true)
    end
    defaults = todict(defaults)
    defaults = Dict{Any, Any}(value(k) => value(v) for (k, v) in pairs(defaults))

    var_to_name = Dict()
    process_variables!(var_to_name, defaults, dvs′)
    process_variables!(var_to_name, defaults, ps′)
    isempty(observed) || collect_var_to_name!(var_to_name, (eq.lhs for eq in observed))

    tgrad = RefValue(EMPTY_TGRAD)
    jac = RefValue{Any}(EMPTY_JAC)
    ctrl_jac = RefValue{Any}(EMPTY_JAC)
    Wfact = RefValue(EMPTY_JAC)
    Wfact_t = RefValue(EMPTY_JAC)
    sysnames = nameof.(systems)
    if length(unique(sysnames)) != length(sysnames)
        throw(ArgumentError("System names must be unique."))
    end
    cont_callbacks = SymbolicContinuousCallbacks(continuous_events)
    disc_callbacks = SymbolicDiscreteCallbacks(discrete_events)
    ODESystem(Threads.atomic_add!(SYSTEM_COUNT, UInt(1)),
        deqs, iv′, dvs′, ps′, tspan, var_to_name, ctrl′, observed, tgrad, jac,
        ctrl_jac, Wfact, Wfact_t, name, systems, defaults, nothing,
        connector_type, preface, cont_callbacks, disc_callbacks,
        metadata, gui_metadata, checks = checks)
end

function ODESystem(eqs, iv = nothing; kwargs...)
    eqs = scalarize(eqs)
    # NOTE: this assumes that the order of algebraic equations doesn't matter
    diffvars = OrderedSet()
    allstates = OrderedSet()
    ps = OrderedSet()
    # reorder equations such that it is in the form of `diffeq, algeeq`
    diffeq = Equation[]
    algeeq = Equation[]
    # initial loop for finding `iv`
    if iv === nothing
        for eq in eqs
            if !(eq.lhs isa Number) # assume eq.lhs is either Differential or Number
                iv = iv_from_nested_derivative(eq.lhs)
                break
            end
        end
    end
    iv = value(iv)
    iv === nothing && throw(ArgumentError("Please pass in independent variables."))
    compressed_eqs = Equation[] # equations that need to be expanded later, like `connect(a, b)`
    for eq in eqs
        eq.lhs isa Union{Symbolic, Number} || (push!(compressed_eqs, eq); continue)
        collect_vars!(allstates, ps, eq.lhs, iv)
        collect_vars!(allstates, ps, eq.rhs, iv)
        if isdiffeq(eq)
            diffvar, _ = var_from_nested_derivative(eq.lhs)
            isequal(iv, iv_from_nested_derivative(eq.lhs)) ||
                throw(ArgumentError("An ODESystem can only have one independent variable."))
            diffvar in diffvars &&
                throw(ArgumentError("The differential variable $diffvar is not unique in the system of equations."))
            push!(diffvars, diffvar)
            push!(diffeq, eq)
        else
            push!(algeeq, eq)
        end
    end
    algevars = setdiff(allstates, diffvars)
    # the orders here are very important!
    return ODESystem(Equation[diffeq; algeeq; compressed_eqs], iv,
        collect(Iterators.flatten((diffvars, algevars))), ps; kwargs...)
end

# NOTE: equality does not check cached Jacobian
function Base.:(==)(sys1::ODESystem, sys2::ODESystem)
    sys1 === sys2 && return true
    iv1 = get_iv(sys1)
    iv2 = get_iv(sys2)
    isequal(iv1, iv2) &&
        isequal(nameof(sys1), nameof(sys2)) &&
        _eq_unordered(get_eqs(sys1), get_eqs(sys2)) &&
        _eq_unordered(get_states(sys1), get_states(sys2)) &&
        _eq_unordered(get_ps(sys1), get_ps(sys2)) &&
        all(s1 == s2 for (s1, s2) in zip(get_systems(sys1), get_systems(sys2)))
end

function flatten(sys::ODESystem, noeqs = false)
    systems = get_systems(sys)
    if isempty(systems)
        return sys
    else
        return ODESystem(noeqs ? Equation[] : equations(sys),
            get_iv(sys),
            states(sys),
            parameters(sys),
            observed = observed(sys),
            continuous_events = continuous_events(sys),
            discrete_events = discrete_events(sys),
            defaults = defaults(sys),
            name = nameof(sys),
            checks = false)
    end
end

ODESystem(eq::Equation, args...; kwargs...) = ODESystem([eq], args...; kwargs...)

"""
$(SIGNATURES)

Build the observed function assuming the observed equations are all explicit,
i.e. there are no cycles.
"""
function build_explicit_observed_function(sys, ts;
    inputs = nothing,
    expression = false,
    output_type = Array,
    checkbounds = true,
    drop_expr = drop_expr,
    throw = true)
    if (isscalar = !(ts isa AbstractVector))
        ts = [ts]
    end
    ts = unwrap.(Symbolics.scalarize(ts))

    vars = Set()
    foreach(Base.Fix1(vars!, vars), ts)
    ivs = independent_variables(sys)
    dep_vars = scalarize(setdiff(vars, ivs))

    obs = observed(sys)

    cs = collect_constants(obs)
    if !isempty(cs) > 0
        cmap = map(x -> x => getdefault(x), cs)
        obs = map(x -> x.lhs ~ substitute(x.rhs, cmap), obs)
    end

    sts = Set(states(sys))
    observed_idx = Dict(x.lhs => i for (i, x) in enumerate(obs))
    param_set = Set(parameters(sys))
    param_set_ns = Set(states(sys, p) for p in parameters(sys))
    namespaced_to_obs = Dict(states(sys, x.lhs) => x.lhs for x in obs)
    namespaced_to_sts = Dict(states(sys, x) => x for x in states(sys))

    # FIXME: This is a rather rough estimate of dependencies. We assume
    # the expression depends on everything before the `maxidx`.
    subs = Dict()
    maxidx = 0
    for s in dep_vars
        if s in param_set || s in param_set_ns
            continue
        end
        idx = get(observed_idx, s, nothing)
        if idx !== nothing
            idx > maxidx && (maxidx = idx)
        else
            s′ = get(namespaced_to_obs, s, nothing)
            if s′ !== nothing
                subs[s] = s′
                s = s′
                idx = get(observed_idx, s, nothing)
            end
            if idx !== nothing
                idx > maxidx && (maxidx = idx)
            elseif !(s in sts)
                s′ = get(namespaced_to_sts, s, nothing)
                if s′ !== nothing
                    subs[s] = s′
                    continue
                end
                if throw
                    Base.throw(ArgumentError("$s is neither an observed nor a state variable."))
                else
                    # TODO: return variables that don't exist in the system.
                    return nothing
                end
            end
            continue
        end
    end
    ts = map(t -> substitute(t, subs), ts)
    obsexprs = []
    for i in 1:maxidx
        eq = obs[i]
        lhs = eq.lhs
        rhs = eq.rhs
        push!(obsexprs, lhs ← rhs)
    end

    pars = parameters(sys)
    if inputs !== nothing
        pars = setdiff(pars, inputs) # Inputs have been converted to parameters by io_preprocessing, remove those from the parameter list
    end
    ps = DestructuredArgs(pars, inbounds = !checkbounds)
    dvs = DestructuredArgs(states(sys), inbounds = !checkbounds)
    if inputs === nothing
        args = [dvs, ps, ivs...]
    else
        ipts = DestructuredArgs(inputs, inbounds = !checkbounds)
        args = [dvs, ipts, ps, ivs...]
    end
    pre = get_postprocess_fbody(sys)

    ex = Func(args, [],
        pre(Let(obsexprs,
            isscalar ? ts[1] : MakeArray(ts, output_type),
            false))) |> toexpr
    expression ? ex : drop_expr(@RuntimeGeneratedFunction(ex))
end

function _eq_unordered(a, b)
    length(a) === length(b) || return false
    n = length(a)
    idxs = Set(1:n)
    for x in a
        idx = findfirst(isequal(x), b)
        idx === nothing && return false
        idx ∈ idxs || return false
        delete!(idxs, idx)
    end
    return true
end

# We have a stand-alone function to convert a `NonlinearSystem` or `ODESystem`
# to an `ODESystem` to connect systems, and we later can reply on
# `structural_simplify` to convert `ODESystem`s to `NonlinearSystem`s.
"""
$(TYPEDSIGNATURES)

Convert a `NonlinearSystem` to an `ODESystem` or converts an `ODESystem` to a
new `ODESystem` with a different independent variable.
"""
function convert_system(::Type{<:ODESystem}, sys, t; name = nameof(sys))
    isempty(observed(sys)) ||
        throw(ArgumentError("`convert_system` cannot handle reduced model (i.e. observed(sys) is non-empty)."))
    t = value(t)
    varmap = Dict()
    sts = states(sys)
    newsts = similar(sts, Any)
    for (i, s) in enumerate(sts)
        if istree(s)
            args = arguments(s)
            length(args) == 1 ||
                throw(InvalidSystemException("Illegal state: $s. The state can have at most one argument like `x(t)`."))
            arg = args[1]
            if isequal(arg, t)
                newsts[i] = s
                continue
            end
            ns = similarterm(s, operation(s), Any[t]; metadata = SymbolicUtils.metadata(s))
            newsts[i] = ns
            varmap[s] = ns
        else
            ns = variable(getname(s); T = FnType)(t)
            newsts[i] = ns
            varmap[s] = ns
        end
    end
    sub = Base.Fix2(substitute, varmap)
    if sys isa AbstractODESystem
        iv = only(independent_variables(sys))
        sub.x[iv] = t # otherwise the Differentials aren't fixed
    end
    neweqs = map(sub, equations(sys))
    defs = Dict(sub(k) => sub(v) for (k, v) in defaults(sys))
    return ODESystem(neweqs, t, newsts, parameters(sys); defaults = defs, name = name,
        checks = false)
end

function Symbolics.substitute(sys::ODESystem, rules::Union{Vector{<:Pair}, Dict})
    rules = todict(map(r -> Symbolics.unwrap(r[1]) => Symbolics.unwrap(r[2]),
        collect(rules)))
    eqs = fast_substitute(equations(sys), rules)
    ODESystem(eqs, get_iv(sys); name = nameof(sys))
end

"""
$(SIGNATURES)

Add accumulation variables for `vars`.
"""
function add_accumulations(sys::ODESystem, vars = states(sys))
    avars = [rename(v, Symbol(:accumulation_, getname(v))) for v in vars]
    add_accumulations(sys, avars .=> vars)
end

"""
$(SIGNATURES)

Add accumulation variables for `vars`. `vars` is a vector of pairs in the form
of

```julia
[cumulative_var1 => x + y, cumulative_var2 => x^2]
```
Then, cumulative variables `cumulative_var1` and `cumulative_var2` that computes
the comulative `x + y` and `x^2` would be added to `sys`.
"""
function add_accumulations(sys::ODESystem, vars::Vector{<:Pair})
    eqs = get_eqs(sys)
    avars = map(first, vars)
    if (ints = intersect(avars, states(sys)); !isempty(ints))
        error("$ints already exist in the system!")
    end
    D = Differential(get_iv(sys))
    @set! sys.eqs = [eqs; Equation[D(a) ~ v[2] for (a, v) in zip(avars, vars)]]
    @set! sys.states = [get_states(sys); avars]
    @set! sys.defaults = merge(get_defaults(sys), Dict(a => 0.0 for a in avars))
end
