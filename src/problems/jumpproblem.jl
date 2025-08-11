@fallback_iip_specialize function JumpProcesses.JumpProblem{iip, spec}(
        sys::System, op, tspan::Union{Tuple, Nothing};
        check_compatibility = true, eval_expression = false, eval_module = @__MODULE__,
        checkbounds = false, cse = true, aggregator = JumpProcesses.NullAggregator(),
        callback = nothing, rng = nothing, kwargs...) where {iip, spec}
    check_complete(sys, JumpProblem)
    check_compatibility && check_compatible_system(JumpProblem, sys)

    has_vrjs = any(x -> x isa VariableRateJump, jumps(sys))
    has_eqs = !isempty(equations(sys))
    has_noise = get_noise_eqs(sys) !== nothing

    if (has_vrjs || has_eqs)
        if has_eqs && has_noise
            prob = SDEProblem{iip, spec}(
                sys, op, tspan; check_compatibility = false,
                build_initializeprob = false, checkbounds, cse, check_length = false,
                kwargs...)
        elseif has_eqs
            prob = ODEProblem{iip, spec}(
                sys, op, tspan; check_compatibility = false,
                build_initializeprob = false, checkbounds, cse, check_length = false,
                kwargs...)
        else
            _, u0,
            p = process_SciMLProblem(EmptySciMLFunction{iip}, sys, op;
                t = tspan === nothing ? nothing : tspan[1],
                check_length = false, build_initializeprob = false, kwargs...)
            observedfun = ObservedFunctionCache(sys; eval_expression, eval_module,
                checkbounds, cse)
            f = (du, u, p, t) -> (du .= 0; nothing)
            df = ODEFunction{true, spec}(f; sys, observed = observedfun)
            prob = ODEProblem{true}(df, u0, tspan, p; kwargs...)
        end
    else
        _f, u0,
        p = process_SciMLProblem(EmptySciMLFunction{iip}, sys, op;
            t = tspan === nothing ? nothing : tspan[1], check_length = false, build_initializeprob = false, cse, kwargs...)
        f = DiffEqBase.DISCRETE_INPLACE_DEFAULT

        observedfun = ObservedFunctionCache(
            sys; eval_expression, eval_module, checkbounds, cse)

        df = DiscreteFunction{true, true}(f; sys = sys, observed = observedfun,
            initialization_data = get(_f.kwargs, :initialization_data, nothing))
        prob = DiscreteProblem(df, u0, tspan, p; kwargs...)
    end

    dvs = unknowns(sys)
    unknowntoid = Dict(value(unknown) => i for (i, unknown) in enumerate(dvs))
    js = jumps(sys)
    invttype = prob.tspan[1] === nothing ? Float64 : typeof(1 / prob.tspan[2])

    # handling parameter substitution and empty param vecs
    p = (prob.p isa DiffEqBase.NullParameters || prob.p === nothing) ? Num[] : prob.p

    majpmapper = JumpSysMajParamMapper(sys, p; jseqs = js, rateconsttype = invttype)
    _majs = Vector{MassActionJump}(filter(x -> x isa MassActionJump, js))
    _crjs = Vector{ConstantRateJump}(filter(x -> x isa ConstantRateJump, js))
    vrjs = Vector{VariableRateJump}(filter(x -> x isa VariableRateJump, js))
    majs = isempty(_majs) ? nothing : assemble_maj(_majs, unknowntoid, majpmapper)
    crjs = ConstantRateJump[assemble_crj(sys, j, unknowntoid; eval_expression, eval_module)
                            for j in _crjs]
    vrjs = VariableRateJump[assemble_vrj(sys, j, unknowntoid; eval_expression, eval_module)
                            for j in vrjs]
    jset = JumpSet(Tuple(vrjs), Tuple(crjs), nothing, majs)

    # dep graphs are only for constant rate jumps
    nonvrjs = ArrayPartition(_majs, _crjs)
    if needs_vartojumps_map(aggregator) || needs_depgraph(aggregator) ||
       (aggregator isa JumpProcesses.NullAggregator)
        jdeps = asgraph(sys; eqs = nonvrjs)
        vdeps = variable_dependencies(sys; eqs = nonvrjs)
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
    cbs = process_events(
        sys; callback, eval_expression, eval_module, op, reset_jumps = true)

    if rng !== nothing
        kwargs = (; kwargs..., rng)
    end
    return JumpProblem(prob, aggregator, jset; dep_graph = jtoj, vartojumps_map = vtoj,
        jumptovars_map = jtov, scale_rates = false, nocopy = true,
        callback = cbs, kwargs...)
end

function check_compatible_system(T::Union{Type{JumpProblem}}, sys::System)
    check_time_dependent(sys, T)
    check_not_dde(sys)
    check_no_cost(sys, T)
    check_no_constraints(sys, T)
    check_has_jumps(sys, T)
    check_is_continuous(sys, T)
end

###################### parameter mapper ###########################
struct JumpSysMajParamMapper{U, V, W}
    paramexprs::U     # the parameter expressions to use for each jump rate constant
    sympars::V        # parameters(sys) from the underlying JumpSystem
    subdict::Any           # mapping from an element of parameters(sys) to its current numerical value
end

function JumpSysMajParamMapper(js::System, p; jseqs = nothing, rateconsttype = Float64)
    eqs = (jseqs === nothing) ? jumps(js) : jseqs
    majs = MassActionJump[x for x in eqs if x isa MassActionJump]
    paramexprs = [maj.scaled_rates for maj in majs]
    psyms = reduce(vcat, reorder_parameters(js); init = [])
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

##### MTK dispatches for Symbolic jumps #####
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

"""
    $(TYPEDSIGNATURES)

Push to `munknowns` the variables modified by jump `jump`. `sts` is the list of unknowns of
the system. Return the modified `munknowns`.
"""
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
