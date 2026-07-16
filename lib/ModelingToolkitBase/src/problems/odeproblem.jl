"""
    generate_ODENLStepData(sys, u0, p, mm, nlstep_compile, nlstep_scc; jac = false)

Generate the NLStep data for implicit ODE solvers. This is a stub that throws an error
if called without ModelingToolkit loaded. The actual implementation is provided by
ModelingToolkit when it is loaded.

When `jac = true`, the analytic Jacobian of the teared inner nonlinear system is
generated symbolically and attached to `nlprob.f.jac`, so NonlinearSolve does not
have to recompute it via AD/FD on every Newton iteration.
"""
function generate_ODENLStepData(sys, u0, p, mm, nlstep_compile, nlstep_scc; jac = false)
    error(
        """
        `nlstep=true` requires ModelingToolkit.jl to be loaded.
        Please add `using ModelingToolkit` to your code before creating an ODEProblem with `nlstep=true`.
        """
    )
end

Base.@nospecializeinfer @fallback_iip_specialize function SciMLBase.ODEFunction{iip, spec}(
        sys::System; @nospecialize(u0 = nothing), @nospecialize(p = nothing), tgrad = false, jac = false,
        t = nothing, eval_expression = false, eval_module = @__MODULE__, sparse = false,
        steady_state = false, checkbounds = false, sparsity = false, @nospecialize(analytic = nothing),
        simplify = false, @nospecialize(initialization_data = nothing), expression = Val{false},
        check_compatibility = true, nlstep = false, nlstep_compile = true, nlstep_scc = false,
        optimize = nothing, kwargs...
    ) where {iip, spec}
    check_complete(sys, ODEFunction)
    check_compatibility && check_compatible_system(ODEFunction, sys)
    if has_any_limited(sys)
        throw(
            ArgumentError(
                "the system contains `limited(...)` nodes that were not stripped; " *
                    "`limited` requires the system to be compiled with `mtkcompile`."
            )
        )
    end

    codegen_opts = GeneratedFunctionOptions(;
        expression, wrap_gfw = Val{true}, eval_expression, eval_module,
        compiler_options = get(kwargs, :compiler_options, CompilerOptions()),
        codegen_function_options = Symbolics.CodegenFunctionOptions(;
            checkbounds, optimize, kwargs...
        )
    )

    f = generate_rhs(sys, codegen_opts)

    if spec === SciMLBase.FunctionWrapperSpecialize && iip
        if u0 === nothing || p === nothing || t === nothing
            error("u0, p, and t must be specified for FunctionWrapperSpecialize on ODEFunction.")
        end
        if expression == Val{true}
            f = :($(SciMLBase.wrapfun_iip)($f, ($u0, $u0, $p, $t)))
        else
            f = SciMLBase.wrapfun_iip(f, (u0, u0, p, t))
        end
    end

    if tgrad
        _tgrad = generate_tgrad(sys, codegen_opts; simplify)
    else
        _tgrad = nothing
    end

    if jac
        _jac = generate_jacobian(sys, codegen_opts; simplify, sparse)
    else
        _jac = nothing
    end

    M = calculate_massmatrix(sys)
    _M = concrete_massmatrix(M; sparse, u0)

    if nlstep
        ode_nlstep = generate_ODENLStepData(sys, u0, p, M, nlstep_compile, nlstep_scc; jac)
    else
        ode_nlstep = nothing
    end

    observedfun = ObservedFunctionCache(
        sys; expression, steady_state, eval_expression, eval_module, checkbounds, optimize
    )

    _W_sparsity = W_sparsity(sys)
    W_prototype = calculate_W_prototype(_W_sparsity; u0, sparse)

    args = (; f)
    kwargs = (;
        sys = sys,
        jac = _jac,
        tgrad = _tgrad,
        mass_matrix = _M,
        jac_prototype = W_prototype,
        observed = observedfun,
        sparsity = sparsity ? _W_sparsity : nothing,
        analytic = analytic,
        initialization_data,
        nlstep_data = ode_nlstep,
    )

    odefn = maybe_codegen_scimlfn(expression, ODEFunction{iip, spec}, args; kwargs...)
    if expression != Val{true} && spec === SciMLBase.AutoSpecialize
        odefn = SciMLBase.widen_bounded_type_params(odefn)
    end
    return odefn
end

Base.@nospecializeinfer @fallback_iip_specialize function SciMLBase.ODEProblem{iip, spec}(
        sys::System, @nospecialize(op), tspan;
        @nospecialize(callback = nothing), check_length = true, eval_expression = false,
        expression = Val{false}, eval_module = @__MODULE__, check_compatibility = true,
        _skip_events = false, kwargs...
    ) where {iip, spec}
    check_complete(sys, ODEProblem)
    check_compatibility && check_compatible_system(ODEProblem, sys)

    _iip = resolve_iip(iip, op)
    if _iip === true
        f, u0, p = process_SciMLProblem(
            ODEFunction{true, spec}, sys, op;
            t = tspan !== nothing ? tspan[1] : tspan, check_length, eval_expression,
            eval_module, expression, check_compatibility, kwargs...
        )
    else
        f, u0, p = process_SciMLProblem(
            ODEFunction{false, spec}, sys, op;
            t = tspan !== nothing ? tspan[1] : tspan, check_length, eval_expression,
            eval_module, expression, check_compatibility, kwargs...
        )
    end

    kwargs = process_kwargs(
        sys; expression, callback, eval_expression, eval_module, op, _skip_events, tspan, kwargs...
    )

    ptype = getmetadata(sys, ProblemTypeCtx, StandardODEProblem())
    args = (; f, u0, tspan, p, ptype)
    maybe_codegen_scimlproblem(expression, ODEProblem{_iip}, args; kwargs...)
end

@fallback_iip_specialize function DiffEqBase.SteadyStateProblem{iip, spec}(
        sys::System, op; check_length = true, check_compatibility = true,
        expression = Val{false}, kwargs...
    ) where {iip, spec}
    check_complete(sys, SteadyStateProblem)
    check_compatibility && check_compatible_system(SteadyStateProblem, sys)

    _iip = resolve_iip(iip, op)
    f, u0,
        p = process_SciMLProblem(
        ODEFunction{_iip}, sys, op;
        steady_state = true, check_length, check_compatibility, expression,
        is_steadystateprob = true, kwargs...
    )

    kwargs = process_kwargs(sys; expression, tspan = (0, Inf), kwargs...)
    args = (; f, u0, p)

    maybe_codegen_scimlproblem(expression, SteadyStateProblem{_iip}, args; kwargs...)
end

function check_compatible_system(
        T::Union{
            Type{ODEFunction}, Type{ODEProblem}, Type{DAEFunction},
            Type{DAEProblem}, Type{SteadyStateProblem},
        },
        sys::System
    )
    check_time_dependent(sys, T)
    check_not_dde(sys)
    check_no_cost(sys, T)
    check_no_constraints(sys, T)
    check_no_jumps(sys, T)
    check_no_noise(sys, T)
    return check_is_continuous(sys, T)
end
