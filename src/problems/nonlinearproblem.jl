@fallback_iip_specialize function SciMLBase.NonlinearFunction{iip, spec}(
        sys::System; u0 = nothing, p = nothing, jac = false,
        eval_expression = false, eval_module = @__MODULE__, sparse = false,
        checkbounds = false, sparsity = false, analytic = nothing,
        simplify = false, cse = true, initialization_data = nothing,
        resid_prototype = nothing, check_compatibility = true, expression = Val{false},
        kwargs...) where {iip, spec}
    check_complete(sys, NonlinearFunction)
    check_compatibility && check_compatible_system(NonlinearFunction, sys)

    f = generate_rhs(sys; expression, wrap_gfw = Val{true},
        eval_expression, eval_module, checkbounds = checkbounds, cse,
        kwargs...)

    if spec === SciMLBase.FunctionWrapperSpecialize && iip
        if u0 === nothing || p === nothing
            error("u0, and p must be specified for FunctionWrapperSpecialize on NonlinearFunction.")
        end
        if expression == Val{true}
            f = :($(SciMLBase.wrapfun_iip)($f, ($u0, $u0, $p)))
        else
            f = SciMLBase.wrapfun_iip(f, (u0, u0, p))
        end
    end

    if jac
        _jac = generate_jacobian(sys; expression,
            wrap_gfw = Val{true}, simplify, sparse, cse, eval_expression, eval_module,
            checkbounds, kwargs...)
    else
        _jac = nothing
    end

    observedfun = ObservedFunctionCache(
        sys; steady_state = false, expression, eval_expression, eval_module, checkbounds,
        cse)

    if sparse
        jac_prototype = similar(calculate_jacobian(sys; sparse), eltype(u0))
    else
        jac_prototype = nothing
    end

    kwargs = (;
        sys = sys,
        jac = _jac,
        observed = observedfun,
        analytic = analytic,
        jac_prototype,
        resid_prototype,
        initialization_data)
    args = (; f)

    return maybe_codegen_scimlfn(expression, NonlinearFunction{iip, spec}, args; kwargs...)
end

@fallback_iip_specialize function SciMLBase.NonlinearProblem{iip, spec}(
        sys::System, op; expression = Val{false},
        check_length = true, check_compatibility = true, kwargs...) where {iip, spec}
    check_complete(sys, NonlinearProblem)
    if is_time_dependent(sys)
        sys = NonlinearSystem(sys)
    end
    check_compatibility && check_compatible_system(NonlinearProblem, sys)

    f, u0,
    p = process_SciMLProblem(NonlinearFunction{iip, spec}, sys, op;
        check_length, check_compatibility, expression, kwargs...)

    kwargs = process_kwargs(sys; kwargs...)
    ptype = getmetadata(sys, ProblemTypeCtx, StandardNonlinearProblem())
    args = (; f, u0, p, ptype)

    return maybe_codegen_scimlproblem(expression, NonlinearProblem{iip}, args; kwargs...)
end

@fallback_iip_specialize function SciMLBase.NonlinearLeastSquaresProblem{iip, spec}(
        sys::System, op; check_length = false,
        check_compatibility = true, expression = Val{false}, kwargs...) where {iip, spec}
    check_complete(sys, NonlinearLeastSquaresProblem)
    check_compatibility && check_compatible_system(NonlinearLeastSquaresProblem, sys)

    f, u0,
    p = process_SciMLProblem(NonlinearFunction{iip}, sys, op;
        check_length, expression, kwargs...)

    kwargs = process_kwargs(sys; kwargs...)
    args = (; f, u0, p)

    return maybe_codegen_scimlproblem(
        expression, NonlinearLeastSquaresProblem{iip}, args; kwargs...)
end

function check_compatible_system(
        T::Union{Type{NonlinearFunction}, Type{NonlinearProblem},
            Type{NonlinearLeastSquaresProblem}}, sys::System)
    check_time_independent(sys, T)
    check_not_dde(sys)
    check_no_cost(sys, T)
    check_no_constraints(sys, T)
    check_no_jumps(sys, T)
    check_no_noise(sys, T)
end

function calculate_resid_prototype(N, u0, p)
    u0ElType = u0 === nothing ? Float64 : eltype(u0)
    if SciMLStructures.isscimlstructure(p)
        u0ElType = promote_type(
            eltype(SciMLStructures.canonicalize(SciMLStructures.Tunable(), p)[1]),
            u0ElType)
    end
    return zeros(u0ElType, N)
end
