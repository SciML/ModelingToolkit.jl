@fallback_iip_specialize function SciMLBase.DDEFunction{iip, spec}(
        sys::System, _d = nothing, u0 = nothing, p = nothing;
        eval_expression = false, eval_module = @__MODULE__, checkbounds = false,
        initialization_data = nothing, cse = true, check_compatibility = true,
        sparse = false, simplify = false, analytic = nothing, kwargs...) where {
        iip, spec}
    check_complete(sys, DDEFunction)
    check_compatibility && check_compatible_system(DDEFunction, sys)

    dvs = unknowns(sys)
    ps = parameters(sys)

    f = generate_rhs(sys, dvs, ps; expression = Val{false},
        eval_expression, eval_module, checkbounds = checkbounds, cse,
        kwargs...)

    if spec === SciMLBase.FunctionWrapperSpecialize && iip
        if u0 === nothing || p === nothing || t === nothing
            error("u0, p, and t must be specified for FunctionWrapperSpecialize on DDEFunction.")
        end
        f = SciMLBase.wrapfun_iip(f, (u0, u0, p, t))
    end

    M = calculate_massmatrix(sys)
    _M = concrete_massmatrix(M; sparse, u0)

    observedfun = ObservedFunctionCache(
        sys; eval_expression, eval_module, checkbounds, cse)

    DDEFunction{iip, spec}(f;
        sys = sys,
        mass_matrix = _M,
        observed = observedfun,
        analytic = analytic,
        initialization_data)
end

@fallback_iip_specialize function SciMLBase.DDEProblem{iip, spec}(
        sys::System, u0map, tspan, parammap = SciMLBase.NullParameters();
        callback = nothing, check_length = true, cse = true, checkbounds = false,
        eval_expression = false, eval_module = @__MODULE__, check_compatibility = true,
        u0_constructor = identity,
        kwargs...) where {iip, spec}
    check_complete(sys, DDEProblem)
    check_compatibility && check_compatible_system(DDEProblem, sys)

    f, u0, p = process_SciMLProblem(DDEFunction{iip, spec}, sys, u0map, parammap;
        t = tspan !== nothing ? tspan[1] : tspan, check_length, cse, checkbounds,
        eval_expression, eval_module, check_compatibility, symbolic_u0 = true, kwargs...)

    h = generate_history(
        sys, u0; expression = Val{false}, cse, eval_expression, eval_module,
        checkbounds)
    u0 = float.(h(p, tspan[1]))
    if u0 !== nothing
        u0 = u0_constructor(u0)
    end

    kwargs = process_kwargs(sys; callback, eval_expression, eval_module, kwargs...)

    # Call `remake` so it runs initialization if it is trivial
    return remake(DDEProblem{iip}(f, u0, h, tspan, p; kwargs...))
end

function check_compatible_system(T::Union{Type{DDEFunction}, Type{DDEProblem}}, sys::System)
    check_time_dependent(sys, T)
    check_is_dde(sys)
    check_no_cost(sys, T)
    check_no_constraints(sys, T)
    check_no_jumps(sys, T)
    check_no_noise(sys, T)
    check_is_continuous(sys, T)
end
