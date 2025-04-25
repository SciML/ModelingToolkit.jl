@fallback_iip_specialize function SciMLBase.SDEFunction{iip, spec}(
        sys::System; u0 = nothing, p = nothing, tgrad = false, jac = false,
        t = nothing, eval_expression = false, eval_module = @__MODULE__, sparse = false,
        steady_state = false, checkbounds = false, sparsity = false, analytic = nothing,
        simplify = false, cse = true, initialization_data = nothing,
        check_compatibility = true, kwargs...) where {iip, spec}
    check_complete(sys, SDEFunction)
    check_compatibility && check_compatible_system(SDEFunction, sys)

    dvs = unknowns(sys)
    ps = parameters(sys)
    f = generate_rhs(sys, dvs, ps; expression = Val{false},
        eval_expression, eval_module, checkbounds = checkbounds, cse,
        kwargs...)
    g = generate_diffusion_function(sys, dvs, ps; expression = Val{false},
        eval_expression, eval_module, checkbounds, cse, kwargs...)

    if spec === SciMLBase.FunctionWrapperSpecialize && iip
        if u0 === nothing || p === nothing || t === nothing
            error("u0, p, and t must be specified for FunctionWrapperSpecialize on ODEFunction.")
        end
        f = SciMLBase.wrapfun_iip(f, (u0, u0, p, t))
    end

    if tgrad
        _tgrad = generate_tgrad(sys, dvs, ps; expression = Val{false},
            simplify, cse, eval_expression, eval_module, checkbounds, kwargs...)
    else
        _tgrad = nothing
    end

    if jac
        _jac = generate_jacobian(sys, dvs, ps; expression = Val{false},
            simplify, sparse, cse, eval_expression, eval_module, checkbounds, kwargs...)
    else
        _jac = nothing
    end

    M = calculate_massmatrix(sys)
    _M = concrete_massmatrix(M; sparse, u0)

    observedfun = ObservedFunctionCache(
        sys; steady_state, eval_expression, eval_module, checkbounds, cse)

    _W_sparsity = W_sparsity(sys)
    W_prototype = calculate_W_prototype(_W_sparsity; u0, sparse)

    SDEFunction{iip, spec}(f, g;
        sys = sys,
        jac = _jac,
        tgrad = _tgrad,
        mass_matrix = _M,
        jac_prototype = W_prototype,
        observed = observedfun,
        sparsity = sparsity ? _W_sparsity : nothing,
        analytic = analytic,
        initialization_data)
end

@fallback_iip_specialize function SciMLBase.SDEProblem{iip, spec}(
        sys::System, u0map, tspan, parammap = SciMLBase.NullParameters();
        callback = nothing, check_length = true, eval_expression = false,
        eval_module = @__MODULE__, check_compatibility = true, sparse = false,
        sparsenoise = sparse, kwargs...) where {iip, spec}
    check_complete(sys, SDEProblem)
    check_compatibility && check_compatible_system(SDEProblem, sys)

    f, u0, p = process_SciMLProblem(SDEFunction{iip, spec}, sys, u0map, parammap;
        t = tspan !== nothing ? tspan[1] : tspan, check_length, eval_expression,
        eval_module, check_compatibility, sparse, kwargs...)

    noise, noise_rate_prototype = calculate_noise_and_rate_prototype(sys, u0; sparsenoise)
    kwargs = process_kwargs(sys; callback, eval_expression, eval_module, kwargs...)
    # Call `remake` so it runs initialization if it is trivial
    return remake(SDEProblem{iip}(f, u0, tspan, p; noise, noise_rate_prototype, kwargs...))
end

function check_compatible_system(T::Union{Type{SDEFunction}, Type{SDEProblem}}, sys::System)
    check_time_dependent(sys, T)
    check_not_dde(sys)
    check_no_cost(sys, T)
    check_no_constraints(sys, T)
    check_no_jumps(sys, T)
    check_has_noise(sys, T)
    check_is_continuous(sys, T)
end

function calculate_noise_and_rate_prototype(sys::System, u0; sparsenoise = false)
    noiseeqs = get_noise_eqs(sys)
    if noiseeqs isa AbstractVector
        # diagonal noise
        noise_rate_prototype = nothing
        noise = nothing
    elseif size(noiseeqs, 2) == 1
        # scalar noise
        noise_rate_prototype = nothing
        noise = WienerProcess(0.0, 0.0, 0.0)
    elseif sparsenoise
        I, J, V = findnz(SparseArrays.sparse(noiseeqs))
        noise_rate_prototype = SparseArrays.sparse(I, J, zero(eltype(u0)))
        noise = nothing
    else
        noise_rate_prototype = zeros(eltype(u0), size(noiseeqs))
        noise = nothing
    end
    return noise, noise_rate_prototype
end
