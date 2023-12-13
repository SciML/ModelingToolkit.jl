const JumpType = Union{VariableRateJump, ConstantRateJump, MassActionJump}

# modifies the expression representing an affect function to
# call reset_aggregated_jumps!(integrator).
# assumes iip
function _reset_aggregator!(expr, integrator)
    if expr isa Symbol
        error("Error, encountered a symbol. This should not happen.")
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

@parameters β γ
@variables t S(t) I(t) R(t)
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
    states::Vector
    """The parameters of the system. Must not contain the independent variable."""
    ps::Vector
    """Array variables."""
    var_to_name::Any
    """Observed states."""
    observed::Vector{Equation}
    """The name of the system."""
    name::Symbol
    """The internal systems. These are required to have unique names."""
    systems::Vector{JumpSystem}
    """
    The default values to use when initial conditions and/or
    parameters are not supplied in `ODEProblem`.
    """
    defaults::Dict
    """
    Type of the system.
    """
    connector_type::Any
    """
    A `Vector{SymbolicDiscreteCallback}` that models events. Symbolic
    analog to `SciMLBase.DiscreteCallback` that executes an affect when a given condition is
    true at the end of an integration step. Note, one must make sure to call
    `reset_aggregated_jumps!(integrator)` if using a custom affect function that changes any
    state value or parameter.
    """
    discrete_events::Vector{SymbolicDiscreteCallback}
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

    function JumpSystem{U}(tag, ap::U, iv, states, ps, var_to_name, observed, name, systems,
            defaults, connector_type, devents,
            metadata = nothing, gui_metadata = nothing,
            complete = false;
            checks::Union{Bool, Int} = true) where {U <: ArrayPartition}
        if checks == true || (checks & CheckComponents) > 0
            check_variables(states, iv)
            check_parameters(ps, iv)
        end
        if checks == true || (checks & CheckUnits) > 0
            all_dimensionless([states; ps; iv]) || check_units(ap, iv)
        end
        new{U}(tag, ap, iv, states, ps, var_to_name, observed, name, systems, defaults,
            connector_type, devents, metadata, gui_metadata, complete)
    end
end

function JumpSystem(eqs, iv, states, ps;
        observed = Equation[],
        systems = JumpSystem[],
        default_u0 = Dict(),
        default_p = Dict(),
        defaults = _merge(Dict(default_u0), Dict(default_p)),
        name = nothing,
        connector_type = nothing,
        checks = true,
        continuous_events = nothing,
        discrete_events = nothing,
        metadata = nothing,
        gui_metadata = nothing,
        kwargs...)
    name === nothing &&
        throw(ArgumentError("The `name` keyword must be provided. Please consider using the `@named` macro"))
    eqs = scalarize(eqs)
    sysnames = nameof.(systems)
    if length(unique(sysnames)) != length(sysnames)
        throw(ArgumentError("System names must be unique."))
    end
    ap = ArrayPartition(MassActionJump[], ConstantRateJump[], VariableRateJump[])
    for eq in eqs
        if eq isa MassActionJump
            push!(ap.x[1], eq)
        elseif eq isa ConstantRateJump
            push!(ap.x[2], eq)
        elseif eq isa VariableRateJump
            push!(ap.x[3], eq)
        else
            error("JumpSystem equations must contain MassActionJumps, ConstantRateJumps, or VariableRateJumps.")
        end
    end
    if !(isempty(default_u0) && isempty(default_p))
        Base.depwarn("`default_u0` and `default_p` are deprecated. Use `defaults` instead.",
            :JumpSystem, force = true)
    end
    defaults = todict(defaults)
    defaults = Dict(value(k) => value(v) for (k, v) in pairs(defaults))

    states, ps = value.(states), value.(ps)
    var_to_name = Dict()
    process_variables!(var_to_name, defaults, states)
    process_variables!(var_to_name, defaults, ps)
    isempty(observed) || collect_var_to_name!(var_to_name, (eq.lhs for eq in observed))
    (continuous_events === nothing) ||
        error("JumpSystems currently only support discrete events.")
    disc_callbacks = SymbolicDiscreteCallbacks(discrete_events)

    JumpSystem{typeof(ap)}(Threads.atomic_add!(SYSTEM_COUNT, UInt(1)),
        ap, value(iv), states, ps, var_to_name, observed, name, systems,
        defaults, connector_type, disc_callbacks, metadata, gui_metadata,
        checks = checks)
end

function generate_rate_function(js::JumpSystem, rate)
    consts = collect_constants(rate)
    if !isempty(consts) # The SymbolicUtils._build_function method of this case doesn't support postprocess_fbody
        csubs = Dict(c => getdefault(c) for c in consts)
        rate = substitute(rate, csubs)
    end
    rf = build_function(rate, states(js), parameters(js),
        get_iv(js),
        expression = Val{true})
end

function generate_affect_function(js::JumpSystem, affect, outputidxs)
    consts = collect_constants(affect)
    if !isempty(consts) # The SymbolicUtils._build_function method of this case doesn't support postprocess_fbody
        csubs = Dict(c => getdefault(c) for c in consts)
        affect = substitute(affect, csubs)
    end
    compile_affect(affect, js, states(js), parameters(js); outputidxs = outputidxs,
        expression = Val{true}, checkvars = false)
end

function assemble_vrj(js, vrj, statetoid)
    rate = drop_expr(@RuntimeGeneratedFunction(generate_rate_function(js, vrj.rate)))
    outputvars = (value(affect.lhs) for affect in vrj.affect!)
    outputidxs = [statetoid[var] for var in outputvars]
    affect = drop_expr(@RuntimeGeneratedFunction(generate_affect_function(js, vrj.affect!,
        outputidxs)))
    VariableRateJump(rate, affect)
end

function assemble_vrj_expr(js, vrj, statetoid)
    rate = generate_rate_function(js, vrj.rate)
    outputvars = (value(affect.lhs) for affect in vrj.affect!)
    outputidxs = ((statetoid[var] for var in outputvars)...,)
    affect = generate_affect_function(js, vrj.affect!, outputidxs)
    quote
        rate = $rate
        affect = $affect
        VariableRateJump(rate, affect)
    end
end

function assemble_crj(js, crj, statetoid)
    rate = drop_expr(@RuntimeGeneratedFunction(generate_rate_function(js, crj.rate)))
    outputvars = (value(affect.lhs) for affect in crj.affect!)
    outputidxs = [statetoid[var] for var in outputvars]
    affect = drop_expr(@RuntimeGeneratedFunction(generate_affect_function(js, crj.affect!,
        outputidxs)))
    ConstantRateJump(rate, affect)
end

function assemble_crj_expr(js, crj, statetoid)
    rate = generate_rate_function(js, crj.rate)
    outputvars = (value(affect.lhs) for affect in crj.affect!)
    outputidxs = ((statetoid[var] for var in outputvars)...,)
    affect = generate_affect_function(js, crj.affect!, outputidxs)
    quote
        rate = $rate
        affect = $affect
        ConstantRateJump(rate, affect)
    end
end

function numericrstoich(mtrs::Vector{Pair{V, W}}, statetoid) where {V, W}
    rs = Vector{Pair{Int, W}}()
    for (wspec, stoich) in mtrs
        spec = value(wspec)
        if !istree(spec) && _iszero(spec)
            push!(rs, 0 => stoich)
        else
            push!(rs, statetoid[spec] => stoich)
        end
    end
    sort!(rs)
    rs
end

function numericnstoich(mtrs::Vector{Pair{V, W}}, statetoid) where {V, W}
    ns = Vector{Pair{Int, W}}()
    for (wspec, stoich) in mtrs
        spec = value(wspec)
        !istree(spec) && _iszero(spec) &&
            error("Net stoichiometry can not have a species labelled 0.")
        push!(ns, statetoid[spec] => stoich)
    end
    sort!(ns)
end

# assemble a numeric MassActionJump from a MT symbolics MassActionJumps
function assemble_maj(majv::Vector{U}, statetoid, pmapper) where {U <: MassActionJump}
    rs = [numericrstoich(maj.reactant_stoch, statetoid) for maj in majv]
    ns = [numericnstoich(maj.net_stoch, statetoid) for maj in majv]
    MassActionJump(rs, ns; param_mapper = pmapper, nocopy = true)
end

"""
```julia
DiffEqBase.DiscreteProblem(sys::JumpSystem, u0map, tspan,
                           parammap = DiffEqBase.NullParameters;
                           use_union = false,
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
dprob = DiscreteProblem(js, u₀map, tspan, parammap)
```
"""
function DiffEqBase.DiscreteProblem(sys::JumpSystem, u0map, tspan::Union{Tuple, Nothing},
        parammap = DiffEqBase.NullParameters();
        checkbounds = false,
        use_union = true,
        kwargs...)
    dvs = states(sys)
    ps = parameters(sys)

    defs = defaults(sys)
    defs = mergedefaults(defs, parammap, ps)
    defs = mergedefaults(defs, u0map, dvs)

    u0 = varmap_to_vars(u0map, dvs; defaults = defs, tofloat = false)
    p = varmap_to_vars(parammap, ps; defaults = defs, tofloat = false, use_union)

    f = DiffEqBase.DISCRETE_INPLACE_DEFAULT

    # just taken from abstractodesystem.jl for ODEFunction def
    obs = observed(sys)
    observedfun = let sys = sys, dict = Dict()
        function generated_observed(obsvar, u, p, t)
            obs = get!(dict, value(obsvar)) do
                build_explicit_observed_function(sys, obsvar; checkbounds = checkbounds)
            end
            obs(u, p, t)
        end
    end

    df = DiscreteFunction{true, true}(f; syms = Symbol.(states(sys)),
        indepsym = Symbol(get_iv(sys)),
        paramsyms = Symbol.(ps), sys = sys,
        observed = observedfun)
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
dprob = DiscreteProblem(js, u₀map, tspan, parammap)
```
"""
struct DiscreteProblemExpr{iip} end

function DiscreteProblemExpr{iip}(sys::JumpSystem, u0map, tspan::Union{Tuple, Nothing},
        parammap = DiffEqBase.NullParameters();
        use_union = false,
        kwargs...) where {iip}
    dvs = states(sys)
    ps = parameters(sys)
    defs = defaults(sys)

    u0 = varmap_to_vars(u0map, dvs; defaults = defs, tofloat = false)
    p = varmap_to_vars(parammap, ps; defaults = defs, tofloat = false, use_union)
    # identity function to make syms works
    quote
        f = DiffEqBase.DISCRETE_INPLACE_DEFAULT
        u0 = $u0
        p = $p
        tspan = $tspan
        df = DiscreteFunction{true, true}(f; syms = $(Symbol.(states(sys))),
            indepsym = $(Symbol(get_iv(sys))),
            paramsyms = $(Symbol.(parameters(sys))))
        DiscreteProblem(df, u0, tspan, p)
    end
end

"""
```julia
DiffEqBase.JumpProblem(js::JumpSystem, prob, aggregator; kwargs...)
```

Generates a JumpProblem from a JumpSystem.

Continuing the example from the [`DiscreteProblem`](@ref) definition:

```julia
jprob = JumpProblem(js, dprob, Direct())
sol = solve(jprob, SSAStepper())
```
"""
function JumpProcesses.JumpProblem(js::JumpSystem, prob, aggregator; callback = nothing,
        kwargs...)
    statetoid = Dict(value(state) => i for (i, state) in enumerate(states(js)))
    eqs = equations(js)
    invttype = prob.tspan[1] === nothing ? Float64 : typeof(1 / prob.tspan[2])

    # handling parameter substitution and empty param vecs
    p = (prob.p isa DiffEqBase.NullParameters || prob.p === nothing) ? Num[] : prob.p

    majpmapper = JumpSysMajParamMapper(js, p; jseqs = eqs, rateconsttype = invttype)
    majs = isempty(eqs.x[1]) ? nothing : assemble_maj(eqs.x[1], statetoid, majpmapper)
    crjs = ConstantRateJump[assemble_crj(js, j, statetoid) for j in eqs.x[2]]
    vrjs = VariableRateJump[assemble_vrj(js, j, statetoid) for j in eqs.x[3]]
    ((prob isa DiscreteProblem) && !isempty(vrjs)) &&
        error("Use continuous problems such as an ODEProblem or a SDEProblem with VariableRateJumps")
    jset = JumpSet(Tuple(vrjs), Tuple(crjs), nothing, majs)

    if needs_vartojumps_map(aggregator) || needs_depgraph(aggregator)
        jdeps = asgraph(js)
        vdeps = variable_dependencies(js)
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
    cbs = process_events(js; callback, postprocess_affect_expr! = _reset_aggregator!)

    JumpProblem(prob, aggregator, jset; dep_graph = jtoj, vartojumps_map = vtoj,
        jumptovars_map = jtov, scale_rates = false, nocopy = true,
        callback = cbs, kwargs...)
end

### Functions to determine which states a jump depends on
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

### Functions to determine which states are modified by a given jump
function modified_states!(mstates, jump::Union{ConstantRateJump, VariableRateJump}, sts)
    for eq in jump.affect!
        st = eq.lhs
        any(isequal(st), sts) && push!(mstates, st)
    end
    mstates
end

function modified_states!(mstates, jump::MassActionJump, sts)
    for (state, stoich) in jump.net_stoch
        any(isequal(state), sts) && push!(mstates, state)
    end
    mstates
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
    psyms = parameters(js)
    paramdict = Dict(value(k) => value(v) for (k, v) in zip(psyms, p))
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

function updateparams!(::JumpSysMajParamMapper{U, V, W},
        params::Nothing) where {U <: AbstractArray, V <: AbstractArray, W}
    nothing
end

# create the initial parameter vector for use in a MassActionJump
function (ratemap::JumpSysMajParamMapper{
        U,
        V,
        W,
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
