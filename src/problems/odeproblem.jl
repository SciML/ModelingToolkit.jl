@fallback_iip_specialize function SciMLBase.ODEFunction{iip, spec}(
        sys::System; u0 = nothing, p = nothing, tgrad = false, jac = false,
        t = nothing, eval_expression = false, eval_module = @__MODULE__, sparse = false,
        steady_state = false, checkbounds = false, sparsity = false, analytic = nothing,
        simplify = false, cse = true, initialization_data = nothing, expression = Val{false},
        check_compatibility = true, nlstep = false, nlstep_compile = true, nlstep_scc = false,
        kwargs...) where {iip, spec}
    check_complete(sys, ODEFunction)
    check_compatibility && check_compatible_system(ODEFunction, sys)

    f = generate_rhs(sys; expression, wrap_gfw = Val{true},
        eval_expression, eval_module, checkbounds = checkbounds, cse,
        kwargs...)

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
        _tgrad = generate_tgrad(
            sys; expression, wrap_gfw = Val{true},
            simplify, cse, eval_expression, eval_module, checkbounds, kwargs...)
    else
        _tgrad = nothing
    end

    if jac
        _jac = generate_jacobian(
            sys; expression, wrap_gfw = Val{true},
            simplify, sparse, cse, eval_expression, eval_module, checkbounds, kwargs...)
    else
        _jac = nothing
    end

    M = calculate_massmatrix(sys)
    _M = concrete_massmatrix(M; sparse, u0)

    if nlstep
        ode_nlstep = generate_ODENLStepData(sys, u0, p, M, nlstep_compile, nlstep_scc)
    else
        ode_nlstep = nothing
    end

    observedfun = ObservedFunctionCache(
        sys; expression, steady_state, eval_expression, eval_module, checkbounds, cse)

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
        nlstep_data = ode_nlstep)

    maybe_codegen_scimlfn(expression, ODEFunction{iip, spec}, args; kwargs...)
end

@fallback_iip_specialize function SciMLBase.ODEProblem{iip, spec}(
        sys::System, op, tspan;
        callback = nothing, check_length = true, eval_expression = false,
        expression = Val{false}, eval_module = @__MODULE__, check_compatibility = true,
        kwargs...) where {iip, spec}
    check_complete(sys, ODEProblem)
    check_compatibility && check_compatible_system(ODEProblem, sys)

    f, u0,
    p = process_SciMLProblem(ODEFunction{iip, spec}, sys, op;
        t = tspan !== nothing ? tspan[1] : tspan, check_length, eval_expression,
        eval_module, expression, check_compatibility, kwargs...)

    kwargs = process_kwargs(
        sys; expression, callback, eval_expression, eval_module, op, kwargs...)

    ptype = getmetadata(sys, ProblemTypeCtx, StandardODEProblem())
    args = (; f, u0, tspan, p, ptype)
    maybe_codegen_scimlproblem(expression, ODEProblem{iip}, args; kwargs...)
end

@fallback_iip_specialize function DiffEqBase.SteadyStateProblem{iip, spec}(
        sys::System, op; check_length = true, check_compatibility = true,
        expression = Val{false}, kwargs...) where {iip, spec}
    check_complete(sys, SteadyStateProblem)
    check_compatibility && check_compatible_system(SteadyStateProblem, sys)

    f, u0,
    p = process_SciMLProblem(ODEFunction{iip}, sys, op;
        steady_state = true, check_length, check_compatibility, expression,
        time_dependent_init = false, kwargs...)

    kwargs = process_kwargs(sys; expression, kwargs...)
    args = (; f, u0, p)

    maybe_codegen_scimlproblem(expression, SteadyStateProblem{iip}, args; kwargs...)
end

function check_compatible_system(
        T::Union{Type{ODEFunction}, Type{ODEProblem}, Type{DAEFunction},
            Type{DAEProblem}, Type{SteadyStateProblem}},
        sys::System)
    check_time_dependent(sys, T)
    check_not_dde(sys)
    check_no_cost(sys, T)
    check_no_constraints(sys, T)
    check_no_jumps(sys, T)
    check_no_noise(sys, T)
    check_is_continuous(sys, T)
end
