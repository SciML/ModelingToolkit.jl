"""
$(TYPEDEF)

A system of ordinary differential equations.

# Fields
$(FIELDS)

# Example

```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters σ ρ β
@variables x(t) y(t) z(t)

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]

@named de = ODESystem(eqs,t,[x,y,z],[σ,ρ,β],tspan=(0, 1000.0))
```
"""
struct ODESystem <: AbstractODESystem
    """
    A tag for the system. If two systems have the same tag, then they are
    structurally identical.
    """
    tag::UInt
    """The ODEs defining the system."""
    eqs::Vector{Equation}
    """Independent variable."""
    iv::BasicSymbolic{Real}
    """
    Dependent (unknown) variables. Must not contain the independent variable.

    N.B.: If `torn_matching !== nothing`, this includes all variables. Actual
    ODE unknowns are determined by the `SelectedState()` entries in `torn_matching`.
    """
    unknowns::Vector
    """Parameter variables. Must not contain the independent variable."""
    ps::Vector
    """Time span."""
    tspan::Union{NTuple{2, Any}, Nothing}
    """Array variables."""
    var_to_name::Any
    """Control parameters (some subset of `ps`)."""
    ctrls::Vector
    """Observed equations."""
    observed::Vector{Equation}
    """System of constraints that must be satisfied by the solution to the system."""
    constraintsystem::Union{Nothing, ConstraintsSystem}
    """A set of expressions defining the costs of the system for optimal control."""
    costs::Vector
    """Takes the cost vector and returns a scalar for optimization."""
    consolidate::Union{Nothing, Function}
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
    Note: this field will not be defined until
    [`generate_factorized_W`](@ref) is called on the system.
    """
    Wfact::RefValue{Matrix{Num}}
    """
    Note: this field will not be defined until
    [`generate_factorized_W`](@ref) is called on the system.
    """
    Wfact_t::RefValue{Matrix{Num}}
    """
    The name of the system.
    """
    name::Symbol
    """
    A description of the system.
    """
    description::String
    """
    The internal systems. These are required to have unique names.
    """
    systems::Vector{ODESystem}
    """
    The default values to use when initial conditions and/or
    parameters are not supplied in `ODEProblem`.
    """
    defaults::Dict
    """
    The guesses to use as the initial conditions for the
    initialization system.
    """
    guesses::Dict
    """
    Tearing result specifying how to solve the system.
    """
    torn_matching::Union{Matching, Nothing}
    """
    The system for performing the initialization.
    """
    initializesystem::Union{Nothing, NonlinearSystem}
    """
    Extra equations to be enforced during the initialization sequence.
    """
    initialization_eqs::Vector{Equation}
    """
    The schedule for the code generation process.
    """
    schedule::Any
    """
    Type of the system.
    """
    connector_type::Any
    """
    Inject assignment statements before the evaluation of the RHS function.
    """
    preface::Any
    """
    A `Vector{SymbolicContinuousCallback}` that model events.
    The integrator will use root finding to guarantee that it steps at each zero crossing.
    """
    continuous_events::Vector{SymbolicContinuousCallback}
    """
    A `Vector{SymbolicDiscreteCallback}` that models events. Symbolic
    analog to `SciMLBase.DiscreteCallback` that executes an affect when a given condition is
    true at the end of an integration step.
    """
    discrete_events::Vector{SymbolicDiscreteCallback}
    """
    Topologically sorted parameter dependency equations, where all symbols are parameters and
    the LHS is a single parameter.
    """
    parameter_dependencies::Vector{Equation}
    """
    Mapping of conditions which should be true throughout the solution process to corresponding error
    messages. These will be added to the equations when calling `debug_system`.
    """
    assertions::Dict{BasicSymbolic, String}
    """
    Metadata for the system, to be used by downstream packages.
    """
    metadata::Any
    """
    Metadata for MTK GUI.
    """
    gui_metadata::Union{Nothing, GUIMetadata}
    """
    A boolean indicating if the given `ODESystem` represents a system of DDEs.
    """
    is_dde::Bool
    """
    A list of points to provide to the solver as tstops. Uses the same syntax as discrete
    events.
    """
    tstops::Vector{Any}
    """
    Cache for intermediate tearing state.
    """
    tearing_state::Any
    """
    Substitutions generated by tearing.
    """
    substitutions::Any
    """
    If false, then `sys.x` no longer performs namespacing.
    """
    namespacing::Bool
    """
    If true, denotes the model will not be modified any further.
    """
    complete::Bool
    """
    Cached data for fast symbolic indexing.
    """
    index_cache::Union{Nothing, IndexCache}
    """
    A list of discrete subsystems.
    """
    discrete_subsystems::Any
    """
    A list of actual unknowns needed to be solved by solvers.
    """
    solved_unknowns::Union{Nothing, Vector{Any}}
    """
    A vector of vectors of indices for the split parameters.
    """
    split_idxs::Union{Nothing, Vector{Vector{Int}}}
    """
    The analysis points removed by transformations, representing connections to be
    ignored. The first element of the tuple analysis points connecting systems and
    the second are ones connecting variables (for the trivial form of `connect`).
    """
    ignored_connections::Union{
        Nothing, Tuple{Vector{IgnoredAnalysisPoint}, Vector{IgnoredAnalysisPoint}}}
    """
    The hierarchical parent system before simplification.
    """
    parent::Any

    function ODESystem(
            tag, deqs, iv, dvs, ps, tspan, var_to_name, ctrls,
            observed, constraints, costs, consolidate, tgrad,
            jac, ctrl_jac, Wfact, Wfact_t, name, description, systems, defaults, guesses,
            torn_matching, initializesystem, initialization_eqs, schedule,
            connector_type, preface, cevents,
            devents, parameter_dependencies, assertions = Dict{BasicSymbolic, String}(),
            metadata = nothing, gui_metadata = nothing, is_dde = false,
            tstops = [], tearing_state = nothing, substitutions = nothing,
            namespacing = true, complete = false, index_cache = nothing,
            discrete_subsystems = nothing, solved_unknowns = nothing,
            split_idxs = nothing, ignored_connections = nothing, parent = nothing;
            checks::Union{Bool, Int} = true)
        if checks == true || (checks & CheckComponents) > 0
            check_independent_variables([iv])
            check_variables(dvs, iv)
            check_parameters(ps, iv)
            check_equations(deqs, iv)
            check_equations(equations(cevents), iv)
            check_subsystems(systems)
        end
        if checks == true || (checks & CheckUnits) > 0
            u = __get_unit_type(dvs, ps, iv)
            check_units(u, deqs)
        end
        new(tag, deqs, iv, dvs, ps, tspan, var_to_name,
            ctrls, observed, constraints, costs, consolidate, tgrad, jac,
            ctrl_jac, Wfact, Wfact_t, name, description, systems, defaults, guesses, torn_matching,
            initializesystem, initialization_eqs, schedule, connector_type, preface,
            cevents, devents, parameter_dependencies, assertions, metadata,
            gui_metadata, is_dde, tstops, tearing_state, substitutions, namespacing,
            complete, index_cache,
            discrete_subsystems, solved_unknowns, split_idxs, ignored_connections, parent)
    end
end

function ODESystem(deqs::AbstractVector{<:Equation}, iv, dvs, ps;
        controls = Num[],
        observed = Equation[],
        constraints = Any[],
        costs = Num[],
        consolidate = nothing,
        systems = ODESystem[],
        tspan = nothing,
        name = nothing,
        description = "",
        default_u0 = Dict(),
        default_p = Dict(),
        defaults = _merge(Dict(default_u0), Dict(default_p)),
        guesses = Dict(),
        initializesystem = nothing,
        initialization_eqs = Equation[],
        schedule = nothing,
        connector_type = nothing,
        preface = nothing,
        continuous_events = nothing,
        discrete_events = nothing,
        parameter_dependencies = Equation[],
        assertions = Dict(),
        checks = true,
        metadata = nothing,
        gui_metadata = nothing,
        is_dde = nothing,
        tstops = [],
        discover_from_metadata = true)
    name === nothing &&
        throw(ArgumentError("The `name` keyword must be provided. Please consider using the `@named` macro"))
    @assert all(control -> any(isequal.(control, ps)), controls) "All controls must also be parameters."

    constraintsystem = nothing
    if !isempty(constraints)
        constraintsystem = process_constraint_system(constraints, dvs, ps, iv)
        for p in parameters(constraintsystem)
            !in(p, Set(ps)) && push!(ps, p)
        end
    end

    if !isempty(costs)
        coststs, costps = process_costs(costs, dvs, ps, iv)
        for p in costps
            !in(p, Set(ps)) && push!(ps, p)
        end
    end
    costs = wrap.(costs)

    iv′ = value(iv)
    ps′ = value.(ps)
    ctrl′ = value.(controls)
    dvs′ = value.(dvs)
    dvs′ = filter(x -> !isdelay(x, iv), dvs′)

    parameter_dependencies, ps′ = process_parameter_dependencies(
        parameter_dependencies, ps′)
    if !(isempty(default_u0) && isempty(default_p))
        Base.depwarn(
            "`default_u0` and `default_p` are deprecated. Use `defaults` instead.",
            :ODESystem, force = true)
    end
    defaults = Dict{Any, Any}(todict(defaults))
    guesses = Dict{Any, Any}(todict(guesses))
    var_to_name = Dict()
    let defaults = discover_from_metadata ? defaults : Dict(),
        guesses = discover_from_metadata ? guesses : Dict()

        process_variables!(var_to_name, defaults, guesses, dvs′)
        process_variables!(var_to_name, defaults, guesses, ps′)
        process_variables!(
            var_to_name, defaults, guesses, [eq.lhs for eq in parameter_dependencies])
        process_variables!(
            var_to_name, defaults, guesses, [eq.rhs for eq in parameter_dependencies])
    end
    defaults = Dict{Any, Any}(value(k) => value(v)
    for (k, v) in pairs(defaults) if v !== nothing)
    guesses = Dict{Any, Any}(value(k) => value(v)
    for (k, v) in pairs(guesses) if v !== nothing)

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

    cont_callbacks, disc_callbacks = create_symbolic_events(
        continuous_events, discrete_events, deqs, iv)
    if is_dde === nothing
        is_dde = _check_if_dde(deqs, iv′, systems)
    end

    if !isempty(systems) && !isnothing(constraintsystem)
        conssystems = ConstraintsSystem[]
        for sys in systems
            cons = get_constraintsystem(sys)
            cons !== nothing && push!(conssystems, cons)
        end
        @set! constraintsystem.systems = conssystems
    end
    costs = wrap.(costs)

    if length(costs) > 1 && isnothing(consolidate)
        error("Must specify a consolidation function for the costs vector.")
    elseif length(costs) == 1 && isnothing(consolidate)
        consolidate = u -> u[1]
    end

    assertions = Dict{BasicSymbolic, Any}(unwrap(k) => v for (k, v) in assertions)

    ODESystem(Threads.atomic_add!(SYSTEM_COUNT, UInt(1)),
        deqs, iv′, dvs′, ps′, tspan, var_to_name, ctrl′, observed,
        constraintsystem, costs, consolidate, tgrad, jac,
        ctrl_jac, Wfact, Wfact_t, name, description, systems,
        defaults, guesses, nothing, initializesystem,
        initialization_eqs, schedule, connector_type, preface, cont_callbacks,
        disc_callbacks, parameter_dependencies, assertions,
        metadata, gui_metadata, is_dde, tstops, checks = checks)
end

function ODESystem(eqs, iv; kwargs...)
    diffvars, allunknowns, ps, eqs = process_equations(eqs, iv)

    for eq in get(kwargs, :parameter_dependencies, Equation[])
        collect_vars!(allunknowns, ps, eq, iv)
    end

    for ssys in get(kwargs, :systems, ODESystem[])
        collect_scoped_vars!(allunknowns, ps, ssys, iv)
    end

    for v in allunknowns
        isdelay(v, iv) || continue
        collect_vars!(allunknowns, ps, arguments(v)[1], iv)
    end

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
            push!(new_ps, p)
        end
    end
    algevars = setdiff(allunknowns, diffvars)

    return ODESystem(eqs, iv, collect(Iterators.flatten((diffvars, algevars))),
        collect(new_ps); kwargs...)
end

# NOTE: equality does not check cached Jacobian
function Base.:(==)(sys1::ODESystem, sys2::ODESystem)
    sys1 === sys2 && return true
    iv1 = get_iv(sys1)
    iv2 = get_iv(sys2)
    isequal(iv1, iv2) &&
        isequal(nameof(sys1), nameof(sys2)) &&
        _eq_unordered(get_eqs(sys1), get_eqs(sys2)) &&
        _eq_unordered(get_unknowns(sys1), get_unknowns(sys2)) &&
        _eq_unordered(get_ps(sys1), get_ps(sys2)) &&
        _eq_unordered(continuous_events(sys1), continuous_events(sys2)) &&
        _eq_unordered(discrete_events(sys1), discrete_events(sys2)) &&
        all(s1 == s2 for (s1, s2) in zip(get_systems(sys1), get_systems(sys2))) &&
        isequal(get_constraintsystem(sys1), get_constraintsystem(sys2)) &&
        _eq_unordered(get_costs(sys1), get_costs(sys2))
end

function flatten(sys::ODESystem, noeqs = false)
    systems = get_systems(sys)
    if isempty(systems)
        return sys
    else
        return ODESystem(noeqs ? Equation[] : equations(sys),
            get_iv(sys),
            unknowns(sys),
            parameters(sys; initial_parameters = true),
            parameter_dependencies = parameter_dependencies(sys),
            guesses = guesses(sys),
            observed = observed(sys),
            continuous_events = continuous_events(sys),
            discrete_events = discrete_events(sys),
            defaults = defaults(sys),
            name = nameof(sys),
            description = description(sys),
            initialization_eqs = initialization_equations(sys),
            assertions = assertions(sys),
            is_dde = is_dde(sys),
            tstops = symbolic_tstops(sys),
            metadata = get_metadata(sys),
            checks = false,
            # without this, any defaults/guesses obtained from metadata that were
            # later removed by the user will be re-added. Right now, we just want to
            # retain `defaults(sys)` as-is.
            discover_from_metadata = false)
    end
end

ODESystem(eq::Equation, args...; kwargs...) = ODESystem([eq], args...; kwargs...)

"""
    build_explicit_observed_function(sys, ts; kwargs...) -> Function(s)

Generates a function that computes the observed value(s) `ts` in the system `sys`, while making the assumption that there are no cycles in the equations.

## Arguments 
- `sys`: The system for which to generate the function
- `ts`: The symbolic observed values whose value should be computed

## Keywords
- `return_inplace = false`: If true and the observed value is a vector, then return both the in place and out of place methods.
- `expression = false`: Generates a Julia `Expr`` computing the observed value if `expression` is true
- `eval_expression = false`: If true and `expression = false`, evaluates the returned function in the module `eval_module`
- `output_type = Array` the type of the array generated by a out-of-place vector-valued function
- `param_only = false` if true, only allow the generated function to access system parameters
- `inputs = Any[]` additional symbolic variables that should be provided to the generated function
- `checkbounds = true` checks bounds if true when destructuring parameters
- `op = Operator` sets the recursion terminator for the walk done by `vars` to identify the variables that appear in `ts`. See the documentation for `vars` for more detail.
- `throw = true` if true, throw an error when generating a function for `ts` that reference variables that do not exist.
- `mkarray`: only used if the output is an array (that is, `!isscalar(ts)`  and `ts` is not a tuple, in which case the result will always be a tuple). Called as `mkarray(ts, output_type)` where `ts` are the expressions to put in the array and `output_type` is the argument of the same name passed to build_explicit_observed_function.
- `cse = true`: Whether to use Common Subexpression Elimination (CSE) to generate a more efficient function.

## Returns

The return value will be either:
* a single function `f_oop` if the input is a scalar or if the input is a Vector but `return_inplace` is false
* the out of place and in-place functions `(f_ip, f_oop)` if `return_inplace` is true and the input is a `Vector`

The function(s) `f_oop` (and potentially `f_ip`) will be:
* `RuntimeGeneratedFunction`s by default,
* A Julia `Expr` if `expression` is true,
* A directly evaluated Julia function in the module `eval_module` if `eval_expression` is true and `expression` is false.

The signatures will be of the form `g(...)` with arguments:

- `output` for in-place functions
- `unknowns` if `param_only` is `false`
- `inputs` if `inputs` is an array of symbolic inputs that should be available in `ts` 
- `p...` unconditionally; note that in the case of `MTKParameters` more than one parameters argument may be present, so it must be splatted
- `t` if the system is time-dependent; for example `NonlinearSystem` will not have `t`

For example, a function `g(op, unknowns, p..., inputs, t)` will be the in-place function generated if `return_inplace` is true, `ts` is a vector, 
an array of inputs `inputs` is given, and `param_only` is false for a time-dependent system.
"""
function build_explicit_observed_function(sys, ts;
        inputs = Any[],
        disturbance_inputs = Any[],
        disturbance_argument = false,
        expression = false,
        eval_expression = false,
        eval_module = @__MODULE__,
        output_type = Array,
        checkbounds = true,
        ps = parameters(sys; initial_parameters = true),
        return_inplace = false,
        param_only = false,
        op = Operator,
        throw = true,
        cse = true,
        mkarray = nothing)
    is_tuple = ts isa Tuple
    if is_tuple
        ts = collect(ts)
        output_type = Tuple
    end

    allsyms = all_symbols(sys)
    if symbolic_type(ts) == NotSymbolic() && ts isa AbstractArray
        ts = map(x -> symbol_to_symbolic(sys, x; allsyms), ts)
    else
        ts = symbol_to_symbolic(sys, ts; allsyms)
    end

    vs = ModelingToolkit.vars(ts; op)
    namespace_subs = Dict()
    ns_map = Dict{Any, Any}(renamespace(sys, obs) => obs for obs in observables(sys))
    for sym in unknowns(sys)
        ns_map[renamespace(sys, sym)] = sym
        if iscall(sym) && operation(sym) === getindex
            ns_map[renamespace(sys, arguments(sym)[1])] = arguments(sym)[1]
        end
    end
    for sym in full_parameters(sys)
        ns_map[renamespace(sys, sym)] = sym
        if iscall(sym) && operation(sym) === getindex
            ns_map[renamespace(sys, arguments(sym)[1])] = arguments(sym)[1]
        end
    end
    allsyms = Set(all_symbols(sys))
    iv = has_iv(sys) ? get_iv(sys) : nothing
    for var in vs
        var = unwrap(var)
        newvar = get(ns_map, var, nothing)
        if newvar !== nothing
            namespace_subs[var] = newvar
            var = newvar
        end
        if throw && !var_in_varlist(var, allsyms, iv)
            Base.throw(ArgumentError("Symbol $var is not present in the system."))
        end
    end
    ts = fast_substitute(ts, namespace_subs)

    obsfilter = if param_only
        if is_split(sys)
            let ic = get_index_cache(sys)
                eq -> !(ContinuousTimeseries() in ic.observed_syms_to_timeseries[eq.lhs])
            end
        else
            Returns(false)
        end
    else
        Returns(true)
    end
    dvs = if param_only
        ()
    else
        (unknowns(sys),)
    end
    if isempty(inputs)
        inputs = ()
    else
        ps = setdiff(ps, inputs) # Inputs have been converted to parameters by io_preprocessing, remove those from the parameter list
        inputs = (inputs,)
    end
    if !isempty(disturbance_inputs)
        # Disturbance inputs may or may not be included as inputs, depending on disturbance_argument
        ps = setdiff(ps, disturbance_inputs)
    end
    if disturbance_argument
        disturbance_inputs = (disturbance_inputs,)
    else
        disturbance_inputs = ()
    end
    ps = reorder_parameters(sys, ps)
    iv = if is_time_dependent(sys)
        (get_iv(sys),)
    else
        ()
    end
    args = (dvs..., inputs..., ps..., iv..., disturbance_inputs...)
    p_start = length(dvs) + length(inputs) + 1
    p_end = length(dvs) + length(inputs) + length(ps)
    fns = build_function_wrapper(
        sys, ts, args...; p_start, p_end, filter_observed = obsfilter,
        output_type, mkarray, try_namespaced = true, expression = Val{true}, cse)
    if fns isa Tuple
        if expression
            return return_inplace ? fns : fns[1]
        end
        oop, iip = eval_or_rgf.(fns; eval_expression, eval_module)
        f = GeneratedFunctionWrapper{(
            p_start + is_dde(sys), length(args) - length(ps) + 1 + is_dde(sys), is_split(sys))}(
            oop, iip)
        return return_inplace ? (f, f) : f
    else
        if expression
            return fns
        end
        f = eval_or_rgf(fns; eval_expression, eval_module)
        f = GeneratedFunctionWrapper{(
            p_start + is_dde(sys), length(args) - length(ps) + 1 + is_dde(sys), is_split(sys))}(
            f, nothing)
        return f
    end
end

function populate_delays(delays::Set, obsexprs, histfn, sys, sym)
    _vars_util = vars(sym)
    for v in _vars_util
        v in delays && continue
        iscall(v) && issym(operation(v)) && (args = arguments(v); length(args) == 1) &&
            iscall(only(args)) || continue

        idx = variable_index(sys, operation(v)(get_iv(sys)))
        idx === nothing && error("Delay term $v is not an unknown in the system")
        push!(delays, v)
        push!(obsexprs, v ← histfn(only(args))[idx])
    end
end

function _eq_unordered(a, b)
    # a and b may be multidimensional
    # e.g. comparing noiseeqs of SDESystem
    a = vec(a)
    b = vec(b)
    length(a) === length(b) || return false
    n = length(a)
    idxs = Set(1:n)
    for x in a
        idx = findfirst(isequal(x), b)
        # loop since there might be multiple identical entries in a/b
        # and while we might have already matched the first there could
        # be a second that is equal to x
        while idx !== nothing && !(idx in idxs)
            idx = findnext(isequal(x), b, idx + 1)
        end
        idx === nothing && return false
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
    sts = unknowns(sys)
    newsts = similar(sts, Any)
    for (i, s) in enumerate(sts)
        if iscall(s)
            args = arguments(s)
            length(args) == 1 ||
                throw(InvalidSystemException("Illegal unknown: $s. The unknown can have at most one argument like `x(t)`."))
            arg = args[1]
            if isequal(arg, t)
                newsts[i] = s
                continue
            end
            ns = maketerm(typeof(s), operation(s), Any[t],
                SymbolicUtils.metadata(s))
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

"""
$(SIGNATURES)

Add accumulation variables for `vars`.
"""
function add_accumulations(sys::ODESystem, vars = unknowns(sys))
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
the cumulative `x + y` and `x^2` would be added to `sys`.
"""
function add_accumulations(sys::ODESystem, vars::Vector{<:Pair})
    eqs = get_eqs(sys)
    avars = map(first, vars)
    if (ints = intersect(avars, unknowns(sys)); !isempty(ints))
        error("$ints already exist in the system!")
    end
    D = Differential(get_iv(sys))
    @set! sys.eqs = [eqs; Equation[D(a) ~ v[2] for (a, v) in zip(avars, vars)]]
    @set! sys.unknowns = [get_unknowns(sys); avars]
    @set! sys.defaults = merge(get_defaults(sys), Dict(a => 0.0 for a in avars))
end

function Base.show(io::IO, mime::MIME"text/plain", sys::ODESystem; hint = true, bold = true)
    # Print general AbstractSystem information
    invoke(Base.show, Tuple{typeof(io), typeof(mime), AbstractSystem},
        io, mime, sys; hint, bold)

    name = nameof(sys)

    # Print initialization equations (unique to ODESystems)
    nini = length(initialization_equations(sys))
    nini > 0 && printstyled(io, "\nInitialization equations ($nini):"; bold)
    nini > 0 && hint && print(io, " see initialization_equations($name)")

    return nothing
end

"""
Build the constraint system for the ODESystem.
"""
function process_constraint_system(
        constraints::Vector, sts, ps, iv; consname = :cons)
    isempty(constraints) && return nothing

    constraintsts = OrderedSet()
    constraintps = OrderedSet()
    for cons in constraints
        collect_vars!(constraintsts, constraintps, cons, iv)
        union!(constraintsts, collect_applied_operators(cons, Differential))
    end

    # Validate the states.
    validate_vars_and_find_ps!(constraintsts, constraintps, sts, iv)

    ConstraintsSystem(
        constraints, collect(constraintsts), collect(constraintps); name = consname)
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
Validate that all the variables in an auxiliary system of the ODESystem (constraint or costs) are 
well-formed states or parameters.
 - Callable/delay variables (e.g. of the form x(0.6) should be unknowns of the system (and have one arg, etc.)
 - Callable/delay parameters should be parameters of the system

Return the set of additional parameters found in the system, e.g. in x(p) ~ 3 then p should be added as a 
parameter of the system.
"""
function validate_vars_and_find_ps!(auxvars, auxps, sysvars, iv)
    sts = Set(sysvars)

    for var in auxvars
        if !iscall(var)
            var ∈ sts ||
                throw(ArgumentError("Time-independent variable $var is not an unknown of the system."))
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

            (isparameter(arg) && !isequal(arg, iv)) && push!(auxps, arg)
        else
            var ∈ sts &&
                @warn "Variable $var has no argument. It will be interpreted as $var($iv), and the constraint will apply to the entire interval."
        end
    end
end

"""
Generate a function that takes a solution object and computes the cost function obtained by coalescing the costs vector.
"""
function generate_cost_function(sys::ODESystem, kwargs...)
    costs = get_costs(sys)
    consolidate = get_consolidate(sys)
    iv = get_iv(sys)

    ps = parameters(sys; initial_parameters = false)
    sts = unknowns(sys)
    np = length(ps)
    ns = length(sts)
    stidxmap = Dict([v => i for (i, v) in enumerate(sts)])
    pidxmap = Dict([v => i for (i, v) in enumerate(ps)])

    @variables sol(..)[1:ns]
    for st in vars(costs)
        x = operation(st)
        t = only(arguments(st))
        idx = stidxmap[x(iv)]

        costs = map(c -> Symbolics.fast_substitute(c, Dict(x(t) => sol(t)[idx])), costs)
    end

    _p = reorder_parameters(sys, ps)
    fs = build_function_wrapper(sys, costs, sol, _p..., t; output_type = Array, kwargs...)
    vc_oop, vc_iip = eval_or_rgf.(fs)

    cost(sol, p, t) = consolidate(vc_oop(sol, p, t))
    return cost
end
