"""
    SymbolicMassActionJump(rate, reactant_stoch, net_stoch; kwargs...)

Construct a `MassActionJump` with `scale_rates = false`, suitable for use in a
`JumpSystem`. The rate expression must already include any combinatorial scaling
(e.g. `k / factorial(n)` for an n-th order homotrimerization reaction).

Returns a `MassActionJump` — this is a convenience constructor, not a new type.
"""
function SymbolicMassActionJump(rate, reactant_stoch, net_stoch; scale_rates = false,
        kwargs...)
    if scale_rates
        throw(ArgumentError(
            "SymbolicMassActionJump requires pre-scaled rate expressions " *
            "(scale_rates = false). scale_rates = true is not supported in " *
            "ModelingToolkitBase."))
    end
    MassActionJump(rate, reactant_stoch, net_stoch; scale_rates = false, kwargs...)
end

@fallback_iip_specialize function JumpProcesses.JumpProblem{iip, spec}(
        sys::System, op, tspan::Union{Tuple, Nothing};
        check_compatibility = true, eval_expression = false, eval_module = @__MODULE__,
        checkbounds = false, cse = true, aggregator = JumpProcesses.NullAggregator(),
        callback = nothing, rng = nothing, save_positions = (true, true), kwargs...
    ) where {iip, spec}
    check_complete(sys, JumpProblem)
    check_compatibility && check_compatible_system(JumpProblem, sys)
    if haskey(kwargs, :tstops)
        throw(ArgumentError(
            "Passing `tstops` directly to `JumpProblem(::System, ...)` is not supported. " *
            "Define tstops on the `System` via the `tstops` keyword instead."))
    end

    has_vrjs = any(x -> x isa VariableRateJump, jumps(sys))
    has_eqs = !isempty(equations(sys))
    has_noise = get_noise_eqs(sys) !== nothing || !isempty(brownians(sys))

    if (has_vrjs || has_eqs)
        if has_eqs && has_noise
            prob = SDEProblem{iip, spec}(
                sys, op, tspan; check_compatibility = false,
                build_initializeprob = false, checkbounds, cse, check_length = false,
                _skip_events = true, _skip_tstops = true, kwargs...
            )
        elseif has_eqs
            prob = ODEProblem{iip, spec}(
                sys, op, tspan; check_compatibility = false,
                build_initializeprob = false, checkbounds, cse, check_length = false,
                _skip_events = true, _skip_tstops = true, kwargs...
            )
        else
            _, u0,
                p = process_SciMLProblem(
                EmptySciMLFunction{iip}, sys, op;
                t = tspan === nothing ? nothing : tspan[1],
                check_length = false, build_initializeprob = false, kwargs...
            )
            observedfun = ObservedFunctionCache(
                sys; eval_expression, eval_module,
                checkbounds, cse
            )
            f = (du, u, p, t) -> (du .= 0; nothing)
            df = ODEFunction{true, spec}(f; sys, observed = observedfun)
            prob = ODEProblem{true}(df, u0, tspan, p; kwargs...)
        end
    else
        _f, u0,
            p = process_SciMLProblem(
            EmptySciMLFunction{iip}, sys, op;
            t = tspan === nothing ? nothing : tspan[1], check_length = false, build_initializeprob = false, cse, kwargs...
        )
        f = DiffEqBase.DISCRETE_INPLACE_DEFAULT

        observedfun = ObservedFunctionCache(
            sys; eval_expression, eval_module, checkbounds, cse
        )

        df = DiscreteFunction{true, true}(
            f; sys = sys, observed = observedfun,
            initialization_data = get(_f.kwargs, :initialization_data, nothing)
        )
        prob = DiscreteProblem(df, u0, tspan, p; kwargs...)
    end

    # Create SymbolicTstops for all paths and forward via JumpProblem kwargs.
    # Inner problems (SDEProblem/ODEProblem) are created with _skip_tstops = true
    # to avoid duplication.
    tstops = SymbolicTstops(sys; eval_expression, eval_module)

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
    crjs = ConstantRateJump[
        assemble_crj(sys, j, unknowntoid; eval_expression, eval_module)
            for j in _crjs
    ]
    vrjs = VariableRateJump[
        assemble_vrj(sys, j, unknowntoid; eval_expression, eval_module)
            for j in vrjs
    ]
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
    # preprocess op to convert Symbol keys to Symbolic using main system before passing
    # to process_events (which may create ImplicitDiscreteProblems for affect subsystems)
    op_processed = operating_point_preprocess(sys, op)
    cbs = process_events(
        sys; callback, eval_expression, eval_module, op = op_processed, reset_jumps = true
    )

    if rng !== nothing
        kwargs = (; kwargs..., rng)
    end
    if tstops !== nothing
        kwargs = (; kwargs..., tstops)
    end
    # MTK requires pre-scaled rate expressions; never ask JumpProcesses to rescale.
    return JumpProblem(
        prob, aggregator, jset; dep_graph = jtoj, vartojumps_map = vtoj,
        jumptovars_map = jtov, scale_rates = false, nocopy = true,
        callback = cbs, save_positions, kwargs...
    )
end

function check_compatible_system(T::Union{Type{JumpProblem}}, sys::System)
    check_time_dependent(sys, T)
    check_not_dde(sys)
    check_no_cost(sys, T)
    check_no_constraints(sys, T)
    check_no_poissonians(sys, T)
    check_has_jumps(sys, T)
    return check_is_continuous(sys, T)
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
    psyms = reduce(vcat, reorder_parameters(js); init = SymbolicT[])
    paramdict = Dict(unwrap(k) => unwrap(v) for (k, v) in zip(psyms, reduce(vcat, p; init = rateconsttype[])))
    return JumpSysMajParamMapper{typeof(paramexprs), typeof(psyms), rateconsttype}(
        paramexprs,
        psyms,
        paramdict
    )
end

function updateparams!(
        ratemap::JumpSysMajParamMapper{U, V, W},
        params
    ) where {U <: AbstractArray, V <: AbstractArray, W}
    for (i, p) in enumerate(params)
        sympar = ratemap.sympars[i]
        ratemap.subdict[sympar] = p
    end
    return nothing
end

function updateparams!(
        ratemap::JumpSysMajParamMapper{U, V, W},
        params::MTKParameters
    ) where {U <: AbstractArray, V <: AbstractArray, W}
    for (i, p) in enumerate(ArrayPartition(params...))
        sympar = ratemap.sympars[i]
        ratemap.subdict[sympar] = p
    end
    return nothing
end

function updateparams!(
        ::JumpSysMajParamMapper{U, V, W},
        params::Nothing
    ) where {U <: AbstractArray, V <: AbstractArray, W}
    return nothing
end

# update a maj with parameter vectors
function (ratemap::JumpSysMajParamMapper{U, V, W})(
        maj::MassActionJump, newparams;
        scale_rates = false,
        kwargs...
    ) where {
        U <: AbstractArray,
        V <: AbstractArray, W,
    }
    updateparams!(ratemap, newparams)
    for i in 1:get_num_majumps(maj)
        maj.scaled_rates[i] = convert(
            W,
            value(
                substitute(
                    ratemap.paramexprs[i],
                    ratemap.subdict; fold = Val(true)
                )
            )
        )
    end
    # No scalerates! call — MTK requires pre-scaled rate expressions.
    # The scale_rates kwarg is accepted but ignored for JumpProcesses API compatibility.
    return nothing
end

# create the initial parameter vector for use in a MassActionJump
function (ratemap::JumpSysMajParamMapper{U, V, W})(params) where {U <: AbstractArray, V <: AbstractArray, W}
    updateparams!(ratemap, params)
    return [
        convert(W, value(substitute(paramexpr, ratemap.subdict; fold = Val(true))))
            for paramexpr in ratemap.paramexprs
    ]
end

##### MTK dispatches for Symbolic jumps #####
eqtype_supports_collect_vars(j::MassActionJump) = true
function collect_vars!(unknowns::OrderedSet{SymbolicT}, parameters::OrderedSet{SymbolicT}, j::MassActionJump, iv::Union{SymbolicT, Nothing}, ::Type{op} = Differential; depth = 0) where {op}
    collect_vars!(unknowns, parameters, j.scaled_rates, iv, op; depth)
    for field in (j.reactant_stoch, j.net_stoch)
        for el in field
            collect_vars!(unknowns, parameters, el, iv, op; depth)
        end
    end
    return nothing
end

eqtype_supports_collect_vars(j::Union{ConstantRateJump, VariableRateJump}) = true
function collect_vars!(
        unknowns::OrderedSet{SymbolicT}, parameters::OrderedSet{SymbolicT}, j::Union{ConstantRateJump, VariableRateJump},
        iv::Union{SymbolicT, Nothing}, ::Type{op} = Differential; depth = 0
    ) where {op}
    collect_vars!(unknowns, parameters, j.rate, iv, op; depth)
    for eq in j.affect!
        (eq isa Equation) && collect_vars!(unknowns, parameters, eq, iv, op; depth)
    end
    return nothing
end

### Functions to determine which unknowns a jump depends on
function SU.search_variables!(dep, jump::Union{ConstantRateJump, VariableRateJump}; kw...)
    jr = unwrap(jump.rate)
    (jr isa SymbolicT) && SU.search_variables!(dep, jr; kw...)
    return dep
end

function SU.search_variables!(dep, jump::MassActionJump; is_atomic = SU.default_is_atomic, kw...)
    sr = unwrap(jump.scaled_rates)
    (sr isa SymbolicT) && SU.search_variables!(dep, sr; is_atomic, kw...)
    for varasop in jump.reactant_stoch
        var = unwrap(varasop[1])
        var isa SymbolicT || continue
        is_atomic(var) && push!(dep, var)
    end
    return dep
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
    return munknowns
end

function modified_unknowns!(munknowns, jump::MassActionJump, sts)
    for (unknown, stoich) in jump.net_stoch
        any(isequal(unknown), sts) && push!(munknowns, unknown)
    end
    return munknowns
end
