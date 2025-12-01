function SciMLBase.IntervalNonlinearFunction(
        sys::System; u0 = nothing, p = nothing, eval_expression = false,
        eval_module = @__MODULE__, expression = Val{false}, checkbounds = false,
        analytic = nothing, cse = true, initialization_data = nothing,
        check_compatibility = true, kwargs...)
    check_complete(sys, IntervalNonlinearFunction)
    check_compatibility && check_compatible_system(IntervalNonlinearFunction, sys)

    f = generate_rhs(sys; expression, wrap_gfw = Val{true},
        scalar = true, eval_expression, eval_module, checkbounds, cse, kwargs...)

    observedfun = ObservedFunctionCache(
        sys; steady_state = false, expression, eval_expression, eval_module, checkbounds,
        cse)

    args = (; f)
    kwargs = (;
        sys = sys,
        observed = observedfun,
        analytic = analytic,
        initialization_data)

    return maybe_codegen_scimlfn(
        expression, IntervalNonlinearFunction{false}, args; kwargs...)
end

function SciMLBase.IntervalNonlinearProblem(
        sys::System, uspan::NTuple{2}, parammap = SciMLBase.NullParameters();
        check_compatibility = true, expression = Val{false}, kwargs...)
    check_complete(sys, IntervalNonlinearProblem)
    check_compatibility && check_compatible_system(IntervalNonlinearProblem, sys)

    u0map = unknowns(sys) .=> uspan[1]
    op = anydict([unknowns(sys)[1] => uspan[1]])
    merge!(op, to_varmap(parammap, parameters(sys)))
    f, u0,
    p = process_SciMLProblem(IntervalNonlinearFunction, sys, op;
        check_compatibility, expression, kwargs...)

    kwargs = process_kwargs(sys; kwargs...)
    args = (; f, uspan, p)
    return maybe_codegen_scimlproblem(expression, IntervalNonlinearProblem, args; kwargs...)
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
