const JUMPPROBLEM_SPECIALIZATION_DOCS = """
Unlike the other problem constructors, `JumpProblem` defaults to
`specialize = SciMLBase.FullSpecialize`. JumpProcesses' aggregators need concrete
parameters: under `SciMLBase.AutoDespecialize`, most of them throw when the jumps are
coupled to an ODE, and the others run much slower. Passing a specialization other than
`SciMLBase.FullSpecialize` is deprecated.
"""

# `tspan` has no default here, unlike every other time-dependent problem constructor:
# `JumpProcesses` already defines `JumpProblem(prob, jumps...)`, so a two-argument method
# with an untyped `op` is ambiguous with it (Aqua reports 14 ambiguities).
"""$(problem_docstring(JumpProcesses.JumpProblem, "inner SciMLFunction", true; init = false, tspan_default = false, extra_body = JUMPPROBLEM_SPECIALIZATION_DOCS))"""
@fallback_iip_specialize SciMLBase.FullSpecialize function JumpProcesses.JumpProblem{iip, spec}(
        sys::System, op, tspan::Union{Tuple, Nothing};
        check_compatibility = true, eval_expression = false, eval_module = @__MODULE__,
        checkbounds = false, aggregator = JumpProcesses.NullAggregator(),
        callback = nothing, rng = nothing, save_positions = (true, true), kwargs...
    ) where {iip, spec}
    check_complete(sys, JumpProblem)
    check_compatibility && check_compatible_system(JumpProblem, sys)
    if spec !== SciMLBase.FullSpecialize
        Base.depwarn(
            "`JumpProblem{iip, $spec}(sys::System, ...)` is deprecated: `JumpProblem` " *
                "for a `System` will only support `SciMLBase.FullSpecialize`, its default. " *
                "Use `JumpProblem(sys, op, tspan)` or `JumpProblem{iip}(sys, op, tspan)`.",
            :JumpProblem
        )
    end
    if haskey(kwargs, :tstops)
        throw(
            ArgumentError(
                "Passing `tstops` directly to `JumpProblem(::System, ...)` is not supported. " *
                    "Define tstops on the `System` via the `tstops` keyword instead."
            )
        )
    end
    # `rng` is kept as a named keyword so it never reaches the wrapped problem.
    if rng !== nothing
        throw(
            ArgumentError(
                "`rng` is no longer a keyword argument for `JumpProblem`. Pass `rng` to " *
                    "`solve` or `init` instead, e.g. `solve(jprob, SSAStepper(); rng = my_rng)`."
            )
        )
    end

    has_vrjs = any(x -> x isa VariableRateJump, jumps(sys))
    has_eqs = !isempty(equations(sys))
    has_noise = get_noise_eqs(sys) !== nothing || !isempty(brownians(sys))

    if (has_vrjs || has_eqs)
        if has_eqs && has_noise
            prob = SDEProblem{iip, spec}(
                sys, op, tspan; check_compatibility = false,
                build_initializeprob = false, checkbounds, check_length = false,
                _skip_events = true, _skip_tstops = true, kwargs...
            )
        elseif has_eqs
            prob = ODEProblem{iip, spec}(
                sys, op, tspan; check_compatibility = false,
                build_initializeprob = false, checkbounds, check_length = false,
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
                checkbounds
            )
            f = (du, u, p, t) -> (du .= 0; nothing)
            df = ODEFunction{true, spec}(f; sys, observed = observedfun)
            prob = ODEProblem{true}(
                df, u0, tspan, p; _solver_kwargs(ODEProblem, filter_kwargs(kwargs))...
            )
        end
    else
        _f, u0,
            p = process_SciMLProblem(
            EmptySciMLFunction{iip}, sys, op;
            t = tspan === nothing ? nothing : tspan[1], check_length = false, build_initializeprob = false, kwargs...
        )
        f = SciMLBase.DISCRETE_INPLACE_DEFAULT

        observedfun = ObservedFunctionCache(
            sys; eval_expression, eval_module, checkbounds
        )

        df = DiscreteFunction{true, true}(
            f; sys = sys, observed = observedfun,
            initialization_data = get(_f.kwargs, :initialization_data, nothing)
        )
        prob = DiscreteProblem(
            df, u0, tspan, p; _solver_kwargs(DiscreteProblem, filter_kwargs(kwargs))...
        )
    end

    # Create SymbolicTstops for all paths and forward via JumpProblem kwargs.
    # Inner problems (SDEProblem/ODEProblem) are created with _skip_tstops = true
    # to avoid duplication.
    tstops = SymbolicTstops(
        sys, GeneratedFunctionOptions(; expression = Val{false}, eval_expression, eval_module)
    )

    dvs = flat_unknowns(sys)
    unknowntoid = Dict(value(unknown) => i for (i, unknown) in enumerate(dvs))
    js = jumps(sys)
    _majs = Vector{MassActionJump}(filter(x -> x isa MassActionJump, js))
    majpmapper = JumpSysMajParamMapper(sys, _majs; eval_expression, eval_module)
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
        sys; callback, eval_expression, eval_module, op = op_processed, reset_jumps = true,
        tspan, mtk_verbosity = _problem_verbosity(kwargs)
    )

    if tstops !== nothing
        kwargs = (; kwargs..., tstops)
    end
    return JumpProblem(
        prob, aggregator, jset; dep_graph = jtoj, vartojumps_map = vtoj,
        jumptovars_map = jtov, nocopy = true,
        callback = cbs, save_positions, _solver_kwargs(JumpProblem, kwargs)...
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
# Fills the working rate constants of the `MassActionJump` built from a system. The rate
# expressions are already scaled, so they are evaluated at the current parameters without
# further scaling. They are compiled once, when the problem is built, so each fill is a
# plain function call; the mapper holds no mutable state, so problems may share it.
struct JumpSysMajParamMapper{F}
    ratefun::F    # ratefun(dest, p) writes each jump's rate constant at parameters p
end

function JumpSysMajParamMapper(
        js::System, majs; eval_expression = false, eval_module = @__MODULE__
    )
    rates = [unwrap(maj.scaled_rates) for maj in majs]
    rps = reorder_parameters(js)
    _, ratefun = build_function_wrapper(
        js, rates, collect(Any, rps),
        BuildFunctionWrapperOptions(;
            p_start = 1, p_end = length(rps),
            codegen_function_options = Symbolics.CodegenFunctionOptions(;
                expression = Val{true}
            )
        )
    )
    ratefun = GeneratedFunctionWrapper{(1, 1, is_split(js))}(
        Val{false}, nothing, ratefun; eval_expression, eval_module
    )
    return JumpSysMajParamMapper(ratefun)
end

# fill `dest` with the rate constants at the parameter values `params`
function (ratemap::JumpSysMajParamMapper)(dest, maj::MassActionJump, params)
    ratemap.ratefun(dest, params)
    return nothing
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
        isin(sts, st) && push!(munknowns, st)
    end
    return munknowns
end

function modified_unknowns!(munknowns, jump::MassActionJump, sts)
    for (unknown, stoich) in jump.net_stoch
        isin(sts, unknown) && push!(munknowns, unknown)
    end
    return munknowns
end

isin(sts, st) = any(isequal(st), sts)
isin(sts::AbstractSet, st) = st in sts
