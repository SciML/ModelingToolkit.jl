@fallback_iip_specialize function SciMLBase.ImplicitDiscreteFunction{iip, spec}(
        sys::System; u0 = nothing, p = nothing, t = nothing, eval_expression = false,
        eval_module = @__MODULE__, expression = Val{false},
        checkbounds = false, analytic = nothing, simplify = false, cse = true,
        initialization_data = nothing, check_compatibility = true, kwargs...) where {
        iip, spec}
    check_complete(sys, ImplicitDiscreteFunction)
    check_compatibility && check_compatible_system(ImplicitDiscreteFunction, sys)

    iv = get_iv(sys)
    dvs = unknowns(sys)
    f = generate_rhs(sys; expression, wrap_gfw = Val{true},
        implicit_dae = true, eval_expression, eval_module, checkbounds = checkbounds, cse,
        override_discrete = true, kwargs...)

    if spec === SciMLBase.FunctionWrapperSpecialize && iip
        if u0 === nothing || p === nothing || t === nothing
            error("u0, p, and t must be specified for FunctionWrapperSpecialize on ImplicitDiscreteFunction.")
        end
        f = SciMLBase.wrapfun_iip(f, (u0, u0, u0, p, t))
    end

    if length(dvs) == length(equations(sys))
        resid_prototype = nothing
    else
        resid_prototype = calculate_resid_prototype(length(equations(sys)), u0, p)
    end

    observedfun = ObservedFunctionCache(
        sys; steady_state = false, expression, eval_expression, eval_module, checkbounds, cse)

    args = (; f)
    kwargs = (;
        sys = sys,
        observed = observedfun,
        analytic = analytic,
        initialization_data,
        resid_prototype)

    return maybe_codegen_scimlfn(
        expression, ImplicitDiscreteFunction{iip, spec}, args; kwargs...)
end

@fallback_iip_specialize function SciMLBase.ImplicitDiscreteProblem{iip, spec}(
        sys::System, op, tspan;
        check_compatibility = true, expression = Val{false}, kwargs...) where {iip, spec}
    check_complete(sys, ImplicitDiscreteProblem)
    check_compatibility && check_compatible_system(ImplicitDiscreteProblem, sys)

    dvs = unknowns(sys)
    op = to_varmap(op, dvs)
    add_toterms!(op; replace = true)
    f, u0,
    p = process_SciMLProblem(
        ImplicitDiscreteFunction{iip, spec}, sys, op;
        t = tspan !== nothing ? tspan[1] : tspan, check_compatibility,
        expression, kwargs...)

    kwargs = process_kwargs(sys; kwargs...)
    args = (; f, u0, tspan, p)
    return maybe_codegen_scimlproblem(
        expression, ImplicitDiscreteProblem{iip}, args; kwargs...)
end

function check_compatible_system(
        T::Union{Type{ImplicitDiscreteFunction}, Type{ImplicitDiscreteProblem}}, sys::System)
    check_time_dependent(sys, T)
    check_not_dde(sys)
    check_no_cost(sys, T)
    check_no_constraints(sys, T)
    check_no_jumps(sys, T)
    check_no_noise(sys, T)
    check_is_discrete(sys, T)
end
