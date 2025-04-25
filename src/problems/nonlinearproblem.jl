@fallback_iip_specialize function SciMLBase.NonlinearFunction{iip, spec}(
        sys::System; u0 = nothing, p = nothing, jac = false,
        eval_expression = false, eval_module = @__MODULE__, sparse = false,
        checkbounds = false, sparsity = false, analytic = nothing,
        simplify = false, cse = true, initialization_data = nothing,
        check_compatibility = true, kwargs...) where {iip, spec}
    check_complete(sys, NonlinearFunction)
    check_compatibility && check_compatible_system(NonlinearFunction, sys)

    dvs = unknowns(sys)
    ps = parameters(sys)
    f = generate_rhs(sys, dvs, ps; expression = Val{false},
        eval_expression, eval_module, checkbounds = checkbounds, cse,
        kwargs...)

    if spec === SciMLBase.FunctionWrapperSpecialize && iip
        if u0 === nothing || p === nothing
            error("u0, and p must be specified for FunctionWrapperSpecialize on NonlinearFunction.")
        end
        f = SciMLBase.wrapfun_iip(f, (u0, u0, p))
    end

    if jac
        _jac = generate_jacobian(sys, dvs, ps; expression = Val{false},
            simplify, sparse, cse, eval_expression, eval_module, checkbounds, kwargs...)
    else
        _jac = nothing
    end

    observedfun = ObservedFunctionCache(
        sys; steady_state = false, eval_expression, eval_module, checkbounds, cse)

    if length(dvs) == length(equations(sys))
        resid_prototype = nothing
    else
        resid_prototype = calculate_resid_prototype(length(equations(sys)), u0, p)
    end

    if sparse
        jac_prototype = similar(calculate_jacobian(sys; sparse), eltype(u0))
    else
        jac_prototype = nothing
    end

    NonlinearFunction{iip, spec}(f;
        sys = sys,
        jac = _jac,
        observed = observedfun,
        analytic = analytic,
        jac_prototype,
        resid_prototype,
        initialization_data)
end

@fallback_iip_specialize function SciMLBase.NonlinearProblem{iip, spec}(
        sys::System, u0map, parammap = SciMLBase.NullParameters();
        check_length = true, check_compatibility = true, kwargs...) where {iip, spec}
    check_complete(sys, NonlinearProblem)
    check_compatibility && check_compatible_system(NonlinearProblem, sys)

    f, u0, p = process_SciMLProblem(NonlinearFunction{iip, spec}, sys, u0map, parammap;
        check_length, check_compatibility, kwargs...)

    kwargs = process_kwargs(sys; kwargs...)
    # Call `remake` so it runs initialization if it is trivial
    return remake(NonlinearProblem{iip}(
        f, u0, p, StandardNonlinearProblem(); kwargs...))
end

@fallback_iip_specialize function SciMLBase.NonlinearLeastSquaresProblem{iip, spec}(
        sys::System, u0map, parammap = DiffEqBase.NullParameters(); check_length = false,
        check_compatibility = true, kwargs...) where {iip, spec}
    check_complete(sys, NonlinearLeastSquaresProblem)
    check_compatibility && check_compatible_system(NonlinearLeastSquaresProblem, sys)

    f, u0, p = process_SciMLProblem(NonlinearFunction{iip}, sys, u0map, parammap;
        check_length, kwargs...)

    kwargs = process_kwargs(sys; kwargs...)
    # Call `remake` so it runs initialization if it is trivial
    return remake(NonlinearLeastSquaresProblem{iip}(f, u0, p; kwargs...))
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
