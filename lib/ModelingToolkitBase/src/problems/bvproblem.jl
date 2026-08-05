# Adapts dynamics generated with an explicit control argument, `f(du, x, u, p, t)`,
# to the stacked decision vector the BVP solvers hand out: each mesh-node vector is
# `[states; controls]`, the controls being extra per-node decision variables with no
# differential equations of their own (`f_prototype` tells the solver how many
# defect rows exist).
struct BVPStackedControlRHS{F} <: Function
    f::F
    nx::Int
end
function (w::BVPStackedControlRHS)(du, xu, p, t)
    w.f(du, @view(xu[1:(w.nx)]), @view(xu[(w.nx + 1):end]), p, t)
    return nothing
end
(w::BVPStackedControlRHS)(xu, p, t) = w.f(@view(xu[1:(w.nx)]), @view(xu[(w.nx + 1):end]), p, t)

"""$(problem_docstring(SciMLBase.BVProblem, ODEFunction, true; init = false, extra_body = BV_EXTRA_BODY))"""
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
    ctrls = inputs(sys)
    op = to_varmap(op, dvs)
    # Optimal-control BVPs (systems carrying a cost) deliberately pin only the
    # operating-point initial conditions — the remaining degrees of freedom are
    # what the cost selects — and their algebraic equations (e.g. lifted output
    # bounds) do not imply a fully determined initial state. Only plain algebraic
    # BVPs skip the guess merge and pin every state.
    is_full_ic = has_alg_eqs(sys) && isempty(get_costs(sys))
    # Systems without algebraic equations should use both fixed values + guesses
    # for initialization.
    _op = is_full_ic ? op : merge(Dict(op), Dict(guesses))

    (; eval_expression, eval_module) = opts.fn_opts.codegen
    checkbounds = opts.fn_opts.codegen.codegen.checkbounds

    if isempty(ctrls)
        inner_opts = maybe_derive_t_from_tspan(opts, tspan)
        inner_opts = setproperties(
            inner_opts; fn_opts = setproperties(inner_opts.fn_opts; check_compatibility = false)
        )
        fode, u0,
            p = process_SciMLProblem(
            ODEFunction{_iip, spec}, sys, _op, inner_opts; options_struct = Val(true), kwargs...
        )
        f_prototype = nothing
    else
        # Like `process_DynamicOptProblem`: a plain `ODEFunction` on a compiled
        # system would freeze the inputs at their operating-point parameter values
        # and silently generate control-free dynamics. Generate with an explicit
        # control argument instead and stack the controls onto each node's decision
        # vector; their initial values come from the processed parameter object.
        E && error("`expression = Val{true}` is not supported for BVProblems with control inputs.")
        fin, u0,
            p = process_SciMLProblem(
            ODEInputFunction{_iip, spec}, sys, _op; guesses, inputs = ctrls,
            t = tspan !== nothing ? tspan[1] : tspan, check_compatibility = false,
            checkbounds, time_dependent_init = false, expression = Val{E}, kwargs...
        )
        fode = BVPStackedControlRHS(fin, length(dvs))
        f_prototype = zeros(eltype(u0), length(dvs))
        # The bare wrapper hides the generated function's mass matrix from
        # BVPFunction's `__has_mass_matrix` probe; forward it explicitly so
        # algebraic (lifted-output) rows keep their zero mass rows.
        mass_matrix = fin.mass_matrix
    end
    fcost = generate_bvp_cost(
        sys,
        GeneratedFunctionOptions(;
            expression = Val{false}, wrap_gfw = Val{false}, eval_expression, eval_module,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; checkbounds)
        )
    )

    stidxmap = Dict([v => i for (i, v) in enumerate(dvs)])
    u0_idxs = is_full_ic ? collect(1:length(dvs)) :
        [stidxmap[k] for (k, v) in op if haskey(stidxmap, k)]
    # The boundary-condition residual is generated from the state-length `u0`;
    # the controls are appended to the problem's decision vector afterwards.
    fbc = generate_boundary_conditions(
        sys, u0, u0_idxs, tspan[1],
        GeneratedFunctionOptions(;
            expression = Val{false}, wrap_gfw = Val{true},
            codegen_function_options = Symbolics.CodegenFunctionOptions(; checkbounds)
        )
    )

    bcresid_prototype = zeros(eltype(u0), length(u0_idxs) + length(constraints(sys)))

    if !isempty(ctrls)
        c0 = collect(eltype(u0), SymbolicIndexingInterface.getp(sys, ctrls)(p))
        u0 = vcat(u0, c0)
    end

    bvpfn = if isempty(ctrls)
        BVPFunction{_iip}(fode, fbc; cost = fcost, f_prototype, bcresid_prototype)
    else
        BVPFunction{_iip}(
            fode, fbc; cost = fcost, f_prototype, bcresid_prototype, mass_matrix
        )
    end

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
