@fallback_iip_specialize function SciMLBase.ImplicitDiscreteFunction{iip, spec}(
        sys::System; u0 = nothing, p = nothing,
        t = nothing, eval_expression = false, eval_module = @__MODULE__,
        checkbounds = false, analytic = nothing, simplify = false, cse = true,
        initialization_data = nothing, check_compatibility = true, kwargs...) where {
        iip, spec}
    check_complete(sys, ImplicitDiscreteFunction)
    check_compatibility && check_compatible_system(ImplicitDiscreteFunction, sys)

    iv = get_iv(sys)
    dvs = unknowns(sys)
    ps = parameters(sys)
    f = generate_rhs(sys, dvs, ps; expression = Val{false}, implicit_dae = true,
        eval_expression, eval_module, checkbounds = checkbounds, cse,
        kwargs...)

    if spec === SciMLBase.FunctionWrapperSpecialize && iip
        if u0 === nothing || p === nothing || t === nothing
            error("u0, p, and t must be specified for FunctionWrapperSpecialize on ImplicitDiscreteFunction.")
        end
        f = SciMLBase.wrapfun_iip(f, (u0, u0, u0, p, t))
    end

    observedfun = ObservedFunctionCache(
        sys; steady_state = false, eval_expression, eval_module, checkbounds, cse)

    ImplicitDiscreteFunction{iip, spec}(f;
        sys = sys,
        observed = observedfun,
        analytic = analytic,
        initialization_data)
end

@fallback_iip_specialize function SciMLBase.ImplicitDiscreteProblem{iip, spec}(
        sys::System, u0map, tspan, parammap = SciMLBase.NullParameters();
        check_compatibility = true, kwargs...) where {iip, spec}
    check_complete(sys, ImplicitDiscreteProblem)
    check_compatibility && check_compatible_system(ImplicitDiscreteProblem, sys)

    dvs = unknowns(sys)
    u0map = to_varmap(u0map, dvs)
    add_toterms!(u0map; replace = true)
    f, u0, p = process_SciMLProblem(
        ImplicitDiscreteFunction{iip, spec}, sys, u0map, parammap;
        t = tspan !== nothing ? tspan[1] : tspan, check_compatibility, kwargs...)

    kwargs = process_kwargs(sys; kwargs...)
    # Call `remake` so it runs initialization if it is trivial
    return remake(ImplicitDiscreteProblem{iip}(f, u0, tspan, p; kwargs...))
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
    check_is_implicit(sys, T, DiscreteProblem)
end
