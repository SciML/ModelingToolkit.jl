function SciMLBase.IntervalNonlinearFunction(
        sys::System, _d = nothing, u0 = nothing, p = nothing;
        eval_expression = false, eval_module = @__MODULE__,
        checkbounds = false, analytic = nothing,
        cse = true, initialization_data = nothing,
        check_compatibility = true, kwargs...)
    check_complete(sys, IntervalNonlinearFunction)
    check_compatibility && check_compatible_system(IntervalNonlinearFunction, sys)

    dvs = unknowns(sys)
    ps = parameters(sys)
    f = generate_rhs(sys, dvs, ps; expression = Val{false}, scalar = true,
        eval_expression, eval_module, checkbounds, cse,
        kwargs...)

    observedfun = ObservedFunctionCache(
        sys; steady_state = false, eval_expression, eval_module, checkbounds, cse)

    IntervalNonlinearFunction{false}(f;
        sys = sys,
        observed = observedfun,
        analytic = analytic,
        initialization_data)
end

function SciMLBase.IntervalNonlinearProblem(
        sys::System, uspan::NTuple{2}, parammap = SciMLBase.NullParameters();
        check_compatibility = true, kwargs...)
    check_complete(sys, IntervalNonlinearProblem)
    check_compatibility && check_compatible_system(IntervalNonlinearProblem, sys)

    u0map = unknowns(sys) .=> uspan[1]
    f, u0, p = process_SciMLProblem(IntervalNonlinearFunction, sys, u0map, parammap;
        check_compatibility, kwargs...)

    kwargs = process_kwargs(sys; kwargs...)
    # Call `remake` so it runs initialization if it is trivial
    return remake(IntervalNonlinearProblem(f, uspan, p; kwargs...))
end

function check_compatible_system(
        T::Union{Type{IntervalNonlinearFunction}, Type{IntervalNonlinearProblem}}, sys::System)
    check_time_independent(sys, T)
    if !isone(length(unknowns(sys)))
        throw(SystemCompatibilityError("""
            `$T` requires a system with a single unknown. Found `$(unknowns(sys))`.
        """))
    end
    if !isone(length(equations(sys)))
        throw(SystemCompatibilityError("""
            `$T` requires a system with a single equation. Found `$(equations(sys))`.
        """))
    end
    check_not_dde(sys)
    check_no_cost(sys, T)
    check_no_constraints(sys, T)
    check_no_jumps(sys, T)
    check_no_noise(sys, T)
end
