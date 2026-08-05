@fallback_iip_specialize function SciMLBase.BVProblem{iip, spec}(
        sys::System, op, tspan;
        check_compatibility = true,
        checkbounds = false, eval_expression = false, eval_module = @__MODULE__,
        expression = Val{false}, guesses = Dict(), callback = nothing,
        kwargs...
    ) where {iip, spec}
    fn_opts = SciMLFunctionOptions(;
        t = tspan !== nothing ? tspan[1] : tspan, eval_expression, eval_module,
        checkbounds, check_compatibility, expression, kwargs...
    )
    opts = SciMLProblemOptions(
        sys;
        fn_opts, guesses, time_dependent_init = false,
        build_initializeprob = supports_initialization(sys),
        circular_dependency_max_cycle_length = length(all_symbols(sys)),
        kwargs...
    )
    return BVProblem{iip, spec}(sys, op, tspan, opts; guesses, callback, kwargs...)
end

"""
    SciMLBase.BVProblem{iip, spec}(sys::System, op, tspan, opts::SciMLProblemOptions; guesses = Dict(), callback = nothing, kwargs...)

Public entry point that builds a `BVProblem` directly from a pre-assembled
[`SciMLProblemOptions`](@ref), bypassing the `kwargs...` wrapper above.

`opts.fn_opts.check_compatibility` governs the `BVProblem`-level compatibility check. The
inner `ODEFunction` build always uses `check_compatibility = false` regardless (it must
skip its own, unrelated ODE-compatibility check) and derives `t` from `tspan[1]` when
`opts.fn_opts.t` isn't already set — both applied via `ConstructionBase.setproperties` on
a copy of `opts`, so a caller who builds `opts` directly (without going through the
keyword wrapper above) gets the same defaulting. `guesses` is kept separate from
`opts.guesses`: the latter is normalized (into a `SymmapT`) for `process_SciMLProblem`'s
own use, while this method needs the raw, unnormalized value to merge into `_op` exactly
as before.
"""
function SciMLBase.BVProblem{iip, spec}(
        sys::System, op, tspan, opts::SciMLProblemOptions{E};
        guesses = Dict(), callback = nothing, kwargs...
    ) where {iip, spec, E}
    check_complete(sys, BVProblem)
    opts.fn_opts.check_compatibility && check_compatible_system(BVProblem, sys)
    isnothing(callback) || error("BVP solvers do not support callbacks.")

    _iip = resolve_iip(iip, op)
    dvs = unknowns(sys)
    op = to_varmap(op, dvs)
    # Systems without algebraic equations should use both fixed values + guesses
    # for initialization.
    _op = has_alg_eqs(sys) ? op : merge(Dict(op), Dict(guesses))

    inner_opts = maybe_derive_t_from_tspan(opts, tspan)
    inner_opts = setproperties(
        inner_opts; fn_opts = setproperties(inner_opts.fn_opts; check_compatibility = false)
    )
    fode, u0,
        p = process_SciMLProblem(
        ODEFunction{_iip, spec}, sys, _op, inner_opts; options_struct = Val(true), kwargs...
    )

    (; eval_expression, eval_module) = opts.fn_opts.codegen
    checkbounds = opts.fn_opts.codegen.codegen.checkbounds
    fcost = generate_bvp_cost(
        sys,
        GeneratedFunctionOptions(;
            expression = Val{false}, wrap_gfw = Val{false}, eval_expression, eval_module,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; checkbounds)
        )
    )

    stidxmap = Dict([v => i for (i, v) in enumerate(dvs)])
    u0_idxs = has_alg_eqs(sys) ? collect(1:length(dvs)) :
        [stidxmap[k] for (k, v) in op if haskey(stidxmap, k)]
    fbc = generate_boundary_conditions(
        sys, u0, u0_idxs, tspan[1],
        GeneratedFunctionOptions(;
            expression = Val{false}, wrap_gfw = Val{true},
            codegen_function_options = Symbolics.CodegenFunctionOptions(; checkbounds)
        )
    )

    n_controls = length(default_codegen_inputs(sys))
    f_prototype = n_controls > 0 ? zeros(eltype(u0), length(dvs) - n_controls) : nothing
    bcresid_prototype = zeros(eltype(u0), length(u0_idxs) + length(constraints(sys)))

    bvpfn = BVPFunction{_iip}(fode, fbc; cost = fcost, f_prototype, bcresid_prototype)

    if (length(constraints(sys)) + length(op) > length(dvs))
        @warn "The BVProblem is overdetermined. The total number of conditions (# constraints + # fixed initial values given by op) exceeds the total number of states. The BVP solvers will default to doing a nonlinear least-squares optimization."
    end

    kwargs = process_kwargs(sys; expression = Val{E}, tspan, kwargs...)
    args = (; bvpfn, u0, tspan, p)

    return maybe_codegen_scimlproblem(Val{E}, BVProblem{_iip}, args; kwargs...)
end

function check_compatible_system(T::Type{BVProblem}, sys::System)
    check_time_dependent(sys, T)
    check_not_dde(sys)
    check_no_jumps(sys, T)
    check_no_noise(sys, T)
    check_is_continuous(sys, T)

    return if !isempty(discrete_events(sys)) || !isempty(continuous_events(sys))
        throw(SystemCompatibilityError("BVP solvers do not support events."))
    end
end
