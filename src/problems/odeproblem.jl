@fallback_iip_specialize function SciMLBase.ODEFunction{iip, spec}(
        sys::System, _d = nothing, u0 = nothing, p = nothing; tgrad = false, jac = false,
        t = nothing, eval_expression = false, eval_module = @__MODULE__, sparse = false,
        steady_state = false, checkbounds = false, sparsity = false, analytic = nothing,
        simplify = false, cse = true, initialization_data = nothing,
        check_compatibility = true, kwargs...) where {iip, spec}
    check_complete(sys, ODEFunction)
    check_compatibility && check_compatible_system(ODEFunction, sys)

    dvs = unknowns(sys)
    ps = parameters(sys)
    f_gen = generate_rhs(sys, dvs, ps; expression = Val{true},
        expression_module = eval_module, checkbounds = checkbounds, cse,
        kwargs...)
    f_oop, f_iip = eval_or_rgf.(f_gen; eval_expression, eval_module)
    f = GeneratedFunctionWrapper{(2, 3, is_split(sys))}(f_oop, f_iip)

    if spec === SciMLBase.FunctionWrapperSpecialize && iip
        if u0 === nothing || p === nothing || t === nothing
            error("u0, p, and t must be specified for FunctionWrapperSpecialize on ODEFunction.")
        end
        f = SciMLBase.wrapfun_iip(f, (u0, u0, p, t))
    end

    if tgrad
        tgrad_gen = generate_tgrad(sys, dvs, ps;
            simplify = simplify,
            expression = Val{true},
            expression_module = eval_module, cse,
            checkbounds = checkbounds, kwargs...)
        tgrad_oop, tgrad_iip = eval_or_rgf.(tgrad_gen; eval_expression, eval_module)
        _tgrad = GeneratedFunctionWrapper{(2, 3, is_split(sys))}(tgrad_oop, tgrad_iip)
    else
        _tgrad = nothing
    end

    if jac
        jac_gen = generate_jacobian(sys, dvs, ps;
            simplify = simplify, sparse = sparse,
            expression = Val{true},
            expression_module = eval_module, cse,
            checkbounds = checkbounds, kwargs...)
        jac_oop, jac_iip = eval_or_rgf.(jac_gen; eval_expression, eval_module)

        _jac = GeneratedFunctionWrapper{(2, 3, is_split(sys))}(jac_oop, jac_iip)
    else
        _jac = nothing
    end

    M = calculate_massmatrix(sys)

    _M = if sparse && !(u0 === nothing || M === I)
        SparseArrays.sparse(M)
    elseif u0 === nothing || M === I
        M
    else
        ArrayInterface.restructure(u0 .* u0', M)
    end

    observedfun = ObservedFunctionCache(
        sys; steady_state, eval_expression, eval_module, checkbounds, cse)

    if sparse
        uElType = u0 === nothing ? Float64 : eltype(u0)
        W_prototype = similar(W_sparsity(sys), uElType)
    else
        W_prototype = nothing
    end

    ODEFunction{iip, spec}(f;
        sys = sys,
        jac = _jac,
        tgrad = _tgrad,
        mass_matrix = _M,
        jac_prototype = W_prototype,
        observed = observedfun,
        sparsity = sparsity ? W_sparsity(sys) : nothing,
        analytic = analytic,
        initialization_data)
end

@fallback_iip_specialize function SciMLBase.ODEProblem{iip, spec}(
        sys::System, u0map, tspan, parammap = SciMLBase.NullParameters();
        callback = nothing, check_length = true, eval_expression = false,
        eval_module = @__MODULE__, check_compatibility = true, kwargs...) where {iip, spec}
    check_complete(sys, ODEProblem)
    check_compatibility && check_compatible_system(ODEProblem, sys)

    f, u0, p = process_SciMLProblem(ODEFunction{iip, spec}, sys, u0map, parammap;
        t = tspan !== nothing ? tspan[1] : tspan,
        check_length, eval_expression, eval_module, kwargs...)
    cbs = process_events(sys; callback, eval_expression, eval_module, kwargs...)

    kwargs = filter_kwargs(kwargs)

    kwargs1 = (;)
    if cbs !== nothing
        kwargs1 = merge(kwargs1, (callback = cbs,))
    end

    tstops = SymbolicTstops(sys; eval_expression, eval_module)
    if tstops !== nothing
        kwargs1 = merge(kwargs1, (; tstops))
    end

    # Call `remake` so it runs initialization if it is trivial
    return remake(ODEProblem{iip}(
        f, u0, tspan, p, StandardODEProblem(); kwargs1..., kwargs...))
end

function check_compatible_system(T::Union{Type{ODEFunction}, Type{ODEProblem}}, sys::System)
    if !is_time_dependent(sys)
        throw(SystemCompatibilityError("""
        `$T` requires a time-dependent system.
        """))
    end

    cost = get_costs(sys)
    if cost isa Vector && !isempty(cost) ||
       cost isa Union{BasicSymbolic, Real} && !_iszero(cost)
        throw(SystemCompatibilityError("""
        `$T` will not optimize solutions of systems that have associated cost \
        functions. Solvers for optimal control problems are forthcoming. In order to \
        bypass this error (e.g. to check the cost of a regular solution), pass \
        `allow_cost = true` into the constructor.
        """))
    end

    if !isempty(constraints(sys))
        throw(SystemCompatibilityError("""
        A system with constraints cannot be used to construct an `$T`. Consider a \
        `BVProblem` instead.
        """))
    end

    if !isempty(jumps(sys))
        throw(SystemCompatibilityError("""
            A system with jumps cannot be used to construct an `$T`. Consider a \
            `JumpProblem` instead.
        """))
    end

    if get_noise_eqs(sys) !== nothing
        throw(SystemCompatibilityError("""
            A system with jumps cannot be used to construct an `$T`. Consider an \
            `SDEProblem` instead.
        """))
    end
end
