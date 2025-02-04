const JumpType = Union{VariableRateJump, ConstantRateJump, MassActionJump}

# modifies the expression representing an affect function to
# call reset_aggregated_jumps!(integrator).
# assumes iip
function _reset_aggregator!(expr, integrator)
    if expr isa Symbol
        error("Error, encountered a symbol. This should not happen.")
    end
    if expr isa LineNumberNode
        return
    end

    if (expr.head == :function)
        _reset_aggregator!(expr.args[end], integrator)
    else
        if expr.args[end] == :nothing
            expr.args[end] = :(reset_aggregated_jumps!($integrator))
            push!(expr.args, :nothing)
        else
            _reset_aggregator!(expr.args[end], integrator)
        end
    end

    nothing
end

"""
$(TYPEDEF)

A system of jump processes.

# Fields
$(FIELDS)

# Example

```julia
using ModelingToolkit, JumpProcesses
using ModelingToolkit: t_nounits as t

@parameters β γ
@variables S(t) I(t) R(t)
rate₁   = β*S*I
affect₁ = [S ~ S - 1, I ~ I + 1]
rate₂   = γ*I
affect₂ = [I ~ I - 1, R ~ R + 1]
j₁      = ConstantRateJump(rate₁,affect₁)
j₂      = ConstantRateJump(rate₂,affect₂)
j₃      = MassActionJump(2*β+γ, [R => 1], [S => 1, R => -1])
@named js      = JumpSystem([j₁,j₂,j₃], t, [S,I,R], [β,γ])
```
"""
struct JumpSystem{U <: ArrayPartition} <: AbstractTimeDependentSystem
    """
    A tag for the system. If two systems have the same tag, then they are
    structurally identical.
    """
    tag::UInt
    """
    The jumps of the system. Allowable types are `ConstantRateJump`,
    `VariableRateJump`, `MassActionJump`.
    """
    eqs::U
    """The independent variable, usually time."""
    iv::Any
    """The dependent variables, representing the state of the system.  Must not contain the independent variable."""
    unknowns::Vector
    """The parameters of the system. Must not contain the independent variable."""
    ps::Vector
    """Array variables."""
    var_to_name::Any
    """Observed variables."""
    observed::Vector{Equation}
    """The name of the system."""
    name::Symbol
    """A description of the system."""
    description::String
    """The internal systems. These are required to have unique names."""
    systems::Vector{JumpSystem}
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
    The system for performing the initialization.
    """
    initializesystem::Union{Nothing, NonlinearSystem}
    """
    Extra equations to be enforced during the initialization sequence.
    """
    initialization_eqs::Vector{Equation}
    """
    Type of the system.
    """
    connector_type::Any
    """
    A `Vector{SymbolicContinuousCallback}` that model events.
    The integrator will use root finding to guarantee that it steps at each zero crossing.
    """
    continuous_events::Vector{SymbolicContinuousCallback}
    """
    A `Vector{SymbolicDiscreteCallback}` that models events. Symbolic
    analog to `SciMLBase.DiscreteCallback` that executes an affect when a given condition is
    true at the end of an integration step. Note, one must make sure to call
    `reset_aggregated_jumps!(integrator)` if using a custom affect function that changes any
    unknown value or parameter.
    """
    discrete_events::Vector{SymbolicDiscreteCallback}
    """
    Topologically sorted parameter dependency equations, where all symbols are parameters and
    the LHS is a single parameter.
    """
    parameter_dependencies::Vector{Equation}
    """
    Metadata for the system, to be used by downstream packages.
    """
    metadata::Any
    """
    Metadata for MTK GUI.
    """
    gui_metadata::Union{Nothing, GUIMetadata}
    """
    If a model `sys` is complete, then `sys.x` no longer performs namespacing.
    """
    complete::Bool
    """
    Cached data for fast symbolic indexing.
    """
    index_cache::Union{Nothing, IndexCache}
    isscheduled::Bool

    function JumpSystem{U}(
            tag, ap::U, iv, unknowns, ps, var_to_name, observed, name, description,
            systems, defaults, guesses, initializesystem, initialization_eqs, connector_type,
            cevents, devents,
            parameter_dependencies, metadata = nothing, gui_metadata = nothing,
            complete = false, index_cache = nothing, isscheduled = false;
            checks::Union{Bool, Int} = true) where {U <: ArrayPartition}
        if checks == true || (checks & CheckComponents) > 0
            check_independent_variables([iv])
            check_variables(unknowns, iv)
            check_parameters(ps, iv)
        end
        if checks == true || (checks & CheckUnits) > 0
            u = __get_unit_type(unknowns, ps, iv)
            check_units(u, ap, iv)
        end
        new{U}(tag, ap, iv, unknowns, ps, var_to_name,
            observed, name, description, systems, defaults, guesses, initializesystem,
            initialization_eqs,
            connector_type, cevents, devents, parameter_dependencies, metadata,
            gui_metadata, complete, index_cache, isscheduled)
    end
end
function JumpSystem(tag, ap, iv, states, ps, var_to_name, args...; kwargs...)
    JumpSystem{typeof(ap)}(tag, ap, iv, states, ps, var_to_name, args...; kwargs...)
end

function JumpSystem(eqs, iv, unknowns, ps;
        observed = Equation[],
        systems = JumpSystem[],
        default_u0 = Dict(),
        default_p = Dict(),
        defaults = _merge(Dict(default_u0), Dict(default_p)),
        guesses = Dict(),
        initializesystem = nothing,
        initialization_eqs = Equation[],
        name = nothing,
        description = "",
        connector_type = nothing,
        checks = true,
        continuous_events = nothing,
        discrete_events = nothing,
        parameter_dependencies = Equation[],
        metadata = nothing,
        gui_metadata = nothing,
        kwargs...)

    # variable processing, similar to ODESystem
    name === nothing &&
        throw(ArgumentError("The `name` keyword must be provided. Please consider using the `@named` macro"))
    iv′ = value(iv)
    us′ = value.(unknowns)
    ps′ = value.(ps)
    parameter_dependencies, ps′ = process_parameter_dependencies(
        parameter_dependencies, ps′)
    if !(isempty(default_u0) && isempty(default_p))
        Base.depwarn(
            "`default_u0` and `default_p` are deprecated. Use `defaults` instead.",
            :JumpSystem, force = true)
    end
    defaults = Dict{Any, Any}(todict(defaults))
    guesses = Dict{Any, Any}(todict(guesses))
    var_to_name = Dict()
    process_variables!(var_to_name, defaults, guesses, us′)
    process_variables!(var_to_name, defaults, guesses, ps′)
    process_variables!(
        var_to_name, defaults, guesses, [eq.lhs for eq in parameter_dependencies])
    process_variables!(
        var_to_name, defaults, guesses, [eq.rhs for eq in parameter_dependencies])
    #! format: off
    defaults = Dict{Any, Any}(value(k) => value(v) for (k, v) in pairs(defaults) if value(v) !== nothing)
    guesses = Dict{Any, Any}(value(k) => value(v) for (k, v) in pairs(guesses) if v !== nothing)
    #! format: on
    isempty(observed) || collect_var_to_name!(var_to_name, (eq.lhs for eq in observed))

    sysnames = nameof.(systems)
    if length(unique(sysnames)) != length(sysnames)
        throw(ArgumentError("System names must be unique."))
    end

    # equation processing
    # this and the treatment of continuous events are the only part 
    # unique to JumpSystems
    eqs = scalarize.(eqs)
    ap = ArrayPartition(
        MassActionJump[], ConstantRateJump[], VariableRateJump[], Equation[])
    for eq in eqs
        if eq isa MassActionJump
            push!(ap.x[1], eq)
        elseif eq isa ConstantRateJump
            push!(ap.x[2], eq)
        elseif eq isa VariableRateJump
            push!(ap.x[3], eq)
        elseif eq isa Equation
            push!(ap.x[4], eq)
        else
            error("JumpSystem equations must contain MassActionJumps, ConstantRateJumps, VariableRateJumps, or Equations.")
        end
    end

    cont_callbacks = SymbolicContinuousCallbacks(continuous_events)
    disc_callbacks = SymbolicDiscreteCallbacks(discrete_events)

    JumpSystem{typeof(ap)}(Threads.atomic_add!(SYSTEM_COUNT, UInt(1)),
        ap, iv′, us′, ps′, var_to_name, observed, name, description, systems,
        defaults, guesses, initializesystem, initialization_eqs, connector_type,
        cont_callbacks, disc_callbacks,
        parameter_dependencies, metadata, gui_metadata, checks = checks)
end

##### MTK dispatches for JumpSystems #####
eqtype_supports_collect_vars(j::MassActionJump) = true
function collect_vars!(unknowns, parameters, j::MassActionJump, iv; depth = 0,
        op = Differential)
    collect_vars!(unknowns, parameters, j.scaled_rates, iv; depth, op)
    for field in (j.reactant_stoch, j.net_stoch)
        for el in field
            collect_vars!(unknowns, parameters, el, iv; depth, op)
        end
    end
    return nothing
end

eqtype_supports_collect_vars(j::Union{ConstantRateJump, VariableRateJump}) = true
function collect_vars!(unknowns, parameters, j::Union{ConstantRateJump, VariableRateJump},
        iv; depth = 0, op = Differential)
    collect_vars!(unknowns, parameters, j.rate, iv; depth, op)
    for eq in j.affect!
        (eq isa Equation) && collect_vars!(unknowns, parameters, eq, iv; depth, op)
    end
    return nothing
end

##########################################

has_massactionjumps(js::JumpSystem) = !isempty(equations(js).x[1])
has_constantratejumps(js::JumpSystem) = !isempty(equations(js).x[2])
has_variableratejumps(js::JumpSystem) = !isempty(equations(js).x[3])
has_equations(js::JumpSystem) = !isempty(equations(js).x[4])

function generate_rate_function(js::JumpSystem, rate)
    consts = collect_constants(rate)
    if !isempty(consts) # The SymbolicUtils._build_function method of this case doesn't support postprocess_fbody
        csubs = Dict(c => getdefault(c) for c in consts)
        rate = substitute(rate, csubs)
    end
    p = reorder_parameters(js, parameters(js))
    rf = build_function_wrapper(js, rate, unknowns(js), p...,
        get_iv(js),
        expression = Val{true})
end

function generate_affect_function(js::JumpSystem, affect, outputidxs)
    consts = collect_constants(affect)
    if !isempty(consts) # The SymbolicUtils._build_function method of this case doesn't support postprocess_fbody
        csubs = Dict(c => getdefault(c) for c in consts)
        affect = substitute(affect, csubs)
    end
    compile_affect(
        affect, nothing, js, unknowns(js), parameters(js); outputidxs = outputidxs,
        expression = Val{true}, checkvars = false)
end

function assemble_vrj(
        js, vrj, unknowntoid; eval_expression = false, eval_module = @__MODULE__)
    rate = eval_or_rgf(generate_rate_function(js, vrj.rate); eval_expression, eval_module)

    outputvars = (value(affect.lhs) for affect in vrj.affect!)
    outputidxs = [unknowntoid[var] for var in outputvars]
    affect = eval_or_rgf(generate_affect_function(js, vrj.affect!, outputidxs);
        eval_expression, eval_module)
    VariableRateJump(rate, affect; save_positions = vrj.save_positions)
end

function assemble_vrj_expr(js, vrj, unknowntoid)
    rate = generate_rate_function(js, vrj.rate)
    outputvars = (value(affect.lhs) for affect in vrj.affect!)
    outputidxs = ((unknowntoid[var] for var in outputvars)...,)
    affect = generate_affect_function(js, vrj.affect!, outputidxs)
    quote
        rate = $rate

        affect = $affect
        VariableRateJump(rate, affect)
    end
end

function assemble_crj(
        js, crj, unknowntoid; eval_expression = false, eval_module = @__MODULE__)
    rate = eval_or_rgf(generate_rate_function(js, crj.rate); eval_expression, eval_module)

    outputvars = (value(affect.lhs) for affect in crj.affect!)
    outputidxs = [unknowntoid[var] for var in outputvars]
    affect = eval_or_rgf(generate_affect_function(js, crj.affect!, outputidxs);
        eval_expression, eval_module)
    ConstantRateJump(rate, affect)
end

function assemble_crj_expr(js, crj, unknowntoid)
    rate = generate_rate_function(js, crj.rate)
    outputvars = (value(affect.lhs) for affect in crj.affect!)
    outputidxs = ((unknowntoid[var] for var in outputvars)...,)
    affect = generate_affect_function(js, crj.affect!, outputidxs)
    quote
        rate = $rate

        affect = $affect
        ConstantRateJump(rate, affect)
    end
end

function numericrstoich(mtrs::Vector{Pair{V, W}}, unknowntoid) where {V, W}
    rs = Vector{Pair{Int, W}}()
    for (wspec, stoich) in mtrs
        spec = value(wspec)
        if !iscall(spec) && _iszero(spec)
            push!(rs, 0 => stoich)
        else
            push!(rs, unknowntoid[spec] => stoich)
        end
    end
    sort!(rs)
    rs
end

function numericnstoich(mtrs::Vector{Pair{V, W}}, unknowntoid) where {V, W}
    ns = Vector{Pair{Int, W}}()
    for (wspec, stoich) in mtrs
        spec = value(wspec)
        !iscall(spec) && _iszero(spec) &&
            error("Net stoichiometry can not have a species labelled 0.")
        push!(ns, unknowntoid[spec] => stoich)
    end
    sort!(ns)
end

# assemble a numeric MassActionJump from a MT symbolics MassActionJumps
function assemble_maj(majv::Vector{U}, unknowntoid, pmapper) where {U <: MassActionJump}
    rs = [numericrstoich(maj.reactant_stoch, unknowntoid) for maj in majv]
    ns = [numericnstoich(maj.net_stoch, unknowntoid) for maj in majv]
    MassActionJump(rs, ns; param_mapper = pmapper, nocopy = true)
end

"""
```julia
DiffEqBase.DiscreteProblem(sys::JumpSystem, u0map, tspan,
                           parammap = DiffEqBase.NullParameters;
                           use_union = true,
                           kwargs...)
```

Generates a blank DiscreteProblem for a pure jump JumpSystem to utilize as
its `prob.prob`. This is used in the case where there are no ODEs
and no SDEs associated with the system.

Continuing the example from the [`JumpSystem`](@ref) definition:

```julia
using DiffEqBase, JumpProcesses
u₀map = [S => 999, I => 1, R => 0]
parammap = [β => 0.1 / 1000, γ => 0.01]
tspan = (0.0, 250.0)
dprob = DiscreteProblem(complete(js), u₀map, tspan, parammap)
```
"""
function DiffEqBase.DiscreteProblem(sys::JumpSystem, u0map, tspan::Union{Tuple, Nothing},
        parammap = DiffEqBase.NullParameters();
        use_union = true,
        eval_expression = false,
        eval_module = @__MODULE__,
        kwargs...)
    if !iscomplete(sys)
        error("A completed `JumpSystem` is required. Call `complete` or `structural_simplify` on the system before creating a `DiscreteProblem`")
    end

    if has_equations(sys) || (!isempty(continuous_events(sys)))
        error("The passed in JumpSystem contains `Equation`s or continuous events, please use a problem type that supports these features, such as ODEProblem.")
    end

    _f, u0, p = process_SciMLProblem(EmptySciMLFunction, sys, u0map, parammap;
        t = tspan === nothing ? nothing : tspan[1], use_union, tofloat = false, check_length = false)
    f = DiffEqBase.DISCRETE_INPLACE_DEFAULT

    observedfun = ObservedFunctionCache(
        sys; eval_expression, eval_module, checkbounds = get(kwargs, :checkbounds, false))

    df = DiscreteFunction{true, true}(f; sys = sys, observed = observedfun,
        initialization_data = get(_f.kwargs, :initialization_data, nothing))
    DiscreteProblem(df, u0, tspan, p; kwargs...)
end

"""
```julia
DiffEqBase.DiscreteProblemExpr(sys::JumpSystem, u0map, tspan,
                               parammap = DiffEqBase.NullParameters; kwargs...)
```

Generates a blank DiscreteProblem for a JumpSystem to utilize as its
solving `prob.prob`. This is used in the case where there are no ODEs
and no SDEs associated with the system.

Continuing the example from the [`JumpSystem`](@ref) definition:

```julia
using DiffEqBase, JumpProcesses
u₀map = [S => 999, I => 1, R => 0]
parammap = [β => 0.1 / 1000, γ => 0.01]
tspan = (0.0, 250.0)
dprob = DiscreteProblem(complete(js), u₀map, tspan, parammap)
```
"""
struct DiscreteProblemExpr{iip} end

function DiscreteProblemExpr{iip}(sys::JumpSystem, u0map, tspan::Union{Tuple, Nothing},
        parammap = DiffEqBase.NullParameters();
        use_union = true,
        kwargs...) where {iip}
    if !iscomplete(sys)
        error("A completed `JumpSystem` is required. Call `complete` or `structural_simplify` on the system before creating a `DiscreteProblemExpr`")
    end

    _, u0, p = process_SciMLProblem(EmptySciMLFunction, sys, u0map, parammap;
        t = tspan === nothing ? nothing : tspan[1], use_union, tofloat = false, check_length = false)
    # identity function to make syms works
    quote
        f = DiffEqBase.DISCRETE_INPLACE_DEFAULT
        u0 = $u0
        p = $p
        sys = $sys
        tspan = $tspan
        df = DiscreteFunction{true, true}(f; sys = sys)
        DiscreteProblem(df, u0, tspan, p)
    end
end

"""
```julia
DiffEqBase.ODEProblem(sys::JumpSystem, u0map, tspan,
                           parammap = DiffEqBase.NullParameters;
                           use_union = true,
                           kwargs...)
```

Generates a blank ODEProblem for a pure jump JumpSystem to utilize as its `prob.prob`. This
is used in the case where there are no ODEs and no SDEs associated with the system but there
are jumps with an explicit time dependency (i.e. `VariableRateJump`s). If no jumps have an
explicit time dependence, i.e. all are `ConstantRateJump`s or `MassActionJump`s then
`DiscreteProblem` should be preferred for performance reasons.

Continuing the example from the [`JumpSystem`](@ref) definition:

```julia
using DiffEqBase, JumpProcesses
u₀map = [S => 999, I => 1, R => 0]
parammap = [β => 0.1 / 1000, γ => 0.01]
tspan = (0.0, 250.0)
oprob = ODEProblem(complete(js), u₀map, tspan, parammap)
```
"""
function DiffEqBase.ODEProblem(sys::JumpSystem, u0map, tspan::Union{Tuple, Nothing},
        parammap = DiffEqBase.NullParameters();
        use_union = false,
        eval_expression = false,
        eval_module = @__MODULE__,
        kwargs...)
    if !iscomplete(sys)
        error("A completed `JumpSystem` is required. Call `complete` or `structural_simplify` on the system before creating a `DiscreteProblem`")
    end

    # forward everything to be an ODESystem but the jumps and discrete events
    if has_equations(sys)
        osys = ODESystem(equations(sys).x[4], get_iv(sys), unknowns(sys), parameters(sys);
            observed = observed(sys), name = nameof(sys), description = description(sys),
            systems = get_systems(sys), defaults = defaults(sys), guesses = guesses(sys),
            parameter_dependencies = parameter_dependencies(sys),
            metadata = get_metadata(sys), gui_metadata = get_gui_metadata(sys))
        osys = complete(osys)
        return ODEProblem(osys, u0map, tspan, parammap; check_length = false, kwargs...)
    else
        _, u0, p = process_SciMLProblem(EmptySciMLFunction, sys, u0map, parammap;
            t = tspan === nothing ? nothing : tspan[1], use_union, tofloat = false,
            check_length = false)
        f = (du, u, p, t) -> (du .= 0; nothing)
        observedfun = ObservedFunctionCache(sys; eval_expression, eval_module,
            checkbounds = get(kwargs, :checkbounds, false))
        df = ODEFunction(f; sys, observed = observedfun)
        return ODEProblem(df, u0, tspan, p; kwargs...)
    end
end

"""
```julia
DiffEqBase.JumpProblem(js::JumpSystem, prob, aggregator; kwargs...)
```

Generates a JumpProblem from a JumpSystem.

Continuing the example from the [`DiscreteProblem`](@ref) definition:

```julia
jprob = JumpProblem(complete(js), dprob, Direct())
sol = solve(jprob, SSAStepper())
```
"""
function JumpProcesses.JumpProblem(js::JumpSystem, prob,
        aggregator = JumpProcesses.NullAggregator(); callback = nothing,
        eval_expression = false, eval_module = @__MODULE__, kwargs...)
    if !iscomplete(js)
        error("A completed `JumpSystem` is required. Call `complete` or `structural_simplify` on the system before creating a `JumpProblem`")
    end
    unknowntoid = Dict(value(unknown) => i for (i, unknown) in enumerate(unknowns(js)))
    eqs = equations(js)
    invttype = prob.tspan[1] === nothing ? Float64 : typeof(1 / prob.tspan[2])

    # handling parameter substitution and empty param vecs
    p = (prob.p isa DiffEqBase.NullParameters || prob.p === nothing) ? Num[] : prob.p

    majpmapper = JumpSysMajParamMapper(js, p; jseqs = eqs, rateconsttype = invttype)
    majs = isempty(eqs.x[1]) ? nothing : assemble_maj(eqs.x[1], unknowntoid, majpmapper)
    crjs = ConstantRateJump[assemble_crj(js, j, unknowntoid; eval_expression, eval_module)
                            for j in eqs.x[2]]
    vrjs = VariableRateJump[assemble_vrj(js, j, unknowntoid; eval_expression, eval_module)
                            for j in eqs.x[3]]
    if prob isa DiscreteProblem
        if (!isempty(vrjs) || has_equations(js) || !isempty(continuous_events(js)))
            error("Use continuous problems such as an ODEProblem or a SDEProblem with VariableRateJumps, coupled differential equations, or continuous events.")
        end
    end
    jset = JumpSet(Tuple(vrjs), Tuple(crjs), nothing, majs)

    # dep graphs are only for constant rate jumps
    nonvrjs = ArrayPartition(eqs.x[1], eqs.x[2])
    if needs_vartojumps_map(aggregator) || needs_depgraph(aggregator) ||
       (aggregator isa JumpProcesses.NullAggregator)
        jdeps = asgraph(js; eqs = nonvrjs)
        vdeps = variable_dependencies(js; eqs = nonvrjs)
        vtoj = jdeps.badjlist
        jtov = vdeps.badjlist
        jtoj = needs_depgraph(aggregator) ? eqeq_dependencies(jdeps, vdeps).fadjlist :
               nothing
    else
        vtoj = nothing
        jtov = nothing
        jtoj = nothing
    end

    # handle events, making sure to reset aggregators in the generated affect functions
    cbs = process_events(js; callback, eval_expression, eval_module,
        postprocess_affect_expr! = _reset_aggregator!)

    JumpProblem(prob, aggregator, jset; dep_graph = jtoj, vartojumps_map = vtoj,
        jumptovars_map = jtov, scale_rates = false, nocopy = true,
        callback = cbs, kwargs...)
end

### Functions to determine which unknowns a jump depends on
function get_variables!(dep, jump::Union{ConstantRateJump, VariableRateJump}, variables)
    jr = value(jump.rate)
    (jr isa Symbolic) && get_variables!(dep, jr, variables)
    dep
end

function get_variables!(dep, jump::MassActionJump, variables)
    sr = value(jump.scaled_rates)
    (sr isa Symbolic) && get_variables!(dep, sr, variables)
    for varasop in jump.reactant_stoch
        any(isequal(varasop[1]), variables) && push!(dep, varasop[1])
    end
    dep
end

### Functions to determine which unknowns are modified by a given jump
function modified_unknowns!(munknowns, jump::Union{ConstantRateJump, VariableRateJump}, sts)
    for eq in jump.affect!
        st = eq.lhs
        any(isequal(st), sts) && push!(munknowns, st)
    end
    munknowns
end

function modified_unknowns!(munknowns, jump::MassActionJump, sts)
    for (unknown, stoich) in jump.net_stoch
        any(isequal(unknown), sts) && push!(munknowns, unknown)
    end
    munknowns
end

###################### parameter mapper ###########################
struct JumpSysMajParamMapper{U, V, W}
    paramexprs::U     # the parameter expressions to use for each jump rate constant
    sympars::V        # parameters(sys) from the underlying JumpSystem
    subdict::Any           # mapping from an element of parameters(sys) to its current numerical value
end

function JumpSysMajParamMapper(js::JumpSystem, p; jseqs = nothing, rateconsttype = Float64)
    eqs = (jseqs === nothing) ? equations(js) : jseqs
    paramexprs = [maj.scaled_rates for maj in eqs.x[1]]
    psyms = reduce(vcat, reorder_parameters(js, parameters(js)); init = [])
    paramdict = Dict(value(k) => value(v) for (k, v) in zip(psyms, vcat(p...)))
    JumpSysMajParamMapper{typeof(paramexprs), typeof(psyms), rateconsttype}(paramexprs,
        psyms,
        paramdict)
end

function updateparams!(ratemap::JumpSysMajParamMapper{U, V, W},
        params) where {U <: AbstractArray, V <: AbstractArray, W}
    for (i, p) in enumerate(params)
        sympar = ratemap.sympars[i]
        ratemap.subdict[sympar] = p
    end
    nothing
end

function updateparams!(ratemap::JumpSysMajParamMapper{U, V, W},
        params::MTKParameters) where {U <: AbstractArray, V <: AbstractArray, W}
    for (i, p) in enumerate(ArrayPartition(params...))
        sympar = ratemap.sympars[i]
        ratemap.subdict[sympar] = p
    end
    nothing
end

function updateparams!(::JumpSysMajParamMapper{U, V, W},
        params::Nothing) where {U <: AbstractArray, V <: AbstractArray, W}
    nothing
end

# create the initial parameter vector for use in a MassActionJump
function (ratemap::JumpSysMajParamMapper{
        U,
        V,
        W
})(params) where {U <: AbstractArray,
        V <: AbstractArray, W}
    updateparams!(ratemap, params)
    [convert(W, value(substitute(paramexpr, ratemap.subdict)))
     for paramexpr in ratemap.paramexprs]
end

# update a maj with parameter vectors
function (ratemap::JumpSysMajParamMapper{U, V, W})(maj::MassActionJump, newparams;
        scale_rates,
        kwargs...) where {U <: AbstractArray,
        V <: AbstractArray, W}
    updateparams!(ratemap, newparams)
    for i in 1:get_num_majumps(maj)
        maj.scaled_rates[i] = convert(W,
            value(substitute(ratemap.paramexprs[i],
                ratemap.subdict)))
    end
    scale_rates && JumpProcesses.scalerates!(maj.scaled_rates, maj.reactant_stoch)
    nothing
end
