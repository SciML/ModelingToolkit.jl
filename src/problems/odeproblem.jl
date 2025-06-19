@fallback_iip_specialize function SciMLBase.ODEFunction{iip, spec}(
        sys::System; u0 = nothing, p = nothing, tgrad = false, jac = false,
        t = nothing, eval_expression = false, eval_module = @__MODULE__, sparse = false,
        steady_state = false, checkbounds = false, sparsity = false, analytic = nothing,
        simplify = false, cse = true, initialization_data = nothing, expression = Val{false},
        check_compatibility = true, kwargs...) where {iip, spec}
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
        initialization_data)

    maybe_codegen_scimlfn(expression, ODEFunction{iip, spec}, args; kwargs...)
end

@fallback_iip_specialize function SciMLBase.ODEProblem{iip, spec}(
        sys::System, op, tspan;
        callback = nothing, check_length = true, eval_expression = false,
        expression = Val{false}, eval_module = @__MODULE__, check_compatibility = true,
        kwargs...) where {iip, spec}
    check_complete(sys, ODEProblem)
    check_compatibility && check_compatible_system(ODEProblem, sys)

    f, u0, p = process_SciMLProblem(ODEFunction{iip, spec}, sys, op;
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

    f, u0, p = process_SciMLProblem(ODEFunction{iip}, sys, op;
        steady_state = true, check_length, check_compatibility, expression,
        time_dependent_init = false, kwargs...)

    kwargs = process_kwargs(sys; expression, kwargs...)
    args = (; f, u0, p)

    maybe_codegen_scimlproblem(expression, SteadyStateProblem{iip}, args; kwargs...)
end

struct SemilinearODEFunction{iip, spec} end

@fallback_iip_specialize function SemilinearODEFunction{iip, specialize}(
        sys::System; u0 = nothing, p = nothing, t = nothing,
        semiquadratic_form = nothing, semiquadratic_jacobian = nothing,
        eval_expression = false, eval_module = @__MODULE__,
        expression = Val{false}, sparse = false, check_compatibility = true,
        jac = false, checkbounds = false, cse = true, initialization_data = nothing,
        analytic = nothing, kwargs...) where {iip, specialize}
    check_complete(sys, SemilinearODEFunction)
    check_compatibility && check_compatible_system(SemilinearODEFunction, sys)

    if semiquadratic_form === nothing
        sys = add_semilinear_parameters(sys)
        semiquadratic_form = calculate_split_form(sys; sparse)
    end

    A, B, x2, C = semiquadratic_form
    M = calculate_massmatrix(sys)
    _M = concrete_massmatrix(M; sparse, u0)

    f1, f2 = generate_semiquadratic_functions(
        sys, A, B, x2, C; expression, wrap_gfw = Val{true},
        eval_expression, eval_module, kwargs...)

    if jac
        semiquadratic_jacobian = @something(semiquadratic_jacobian,
            calculate_semiquadratic_jacobian(sys, B, x2, C; sparse, massmatrix = _M))
        f1jac, x2jac, Cjac = semiquadratic_jacobian
        _jac = generate_semiquadratic_jacobian(
            sys, B, x2, C, f1jac, x2jac, Cjac; sparse, expression,
            wrap_gfw = Val{true}, eval_expression, eval_module, kwargs...)
        _W_sparsity = f1jac
        W_prototype = calculate_W_prototype(_W_sparsity; u0, sparse)
    else
        _jac = nothing
        W_prototype = nothing
    end

    observedfun = ObservedFunctionCache(
        sys; expression, steady_state = false, eval_expression, eval_module, checkbounds, cse)

    f1_args = (; f1)
    f1_kwargs = (; jac = _jac)
    f1 = maybe_codegen_scimlfn(
        expression, ODEFunction{iip, specialize}, f1_args; f1_kwargs...)
    args = (; f1, f2)

    kwargs = (;
        sys = sys,
        jac = _jac,
        mass_matrix = _M,
        jac_prototype = W_prototype,
        observed = observedfun,
        analytic,
        initialization_data)
    kwargs = (; sys, observed = observedfun, mass_matrix = _M)

    return maybe_codegen_scimlfn(
        expression, SplitFunction{iip, specialize}, args; kwargs...)
end

struct SemilinearODEProblem{iip, spec} end

@fallback_iip_specialize function SemilinearODEProblem{iip, spec}(
        sys::System, op, tspan; check_compatibility = true,
        u0_eltype = nothing, expression = Val{false}, callback = nothing,
        jac = false, sparse = false, kwargs...) where {iip, spec}
    check_complete(sys, SemilinearODEProblem)
    check_compatibility && check_compatible_system(SemilinearODEProblem, sys)

    A, B, x2, C = semiquadratic_form = calculate_split_form(sys)

    semiquadratic_jacobian = nothing
    if jac
        f1jac, x2jac, Cjac = semiquadratic_jacobian = calculate_semiquadratic_jacobian(
            sys, B, x2, C; sparse)
    end

    sys = add_semilinear_parameters(sys)
    linear_matrix_param = unwrap(getproperty(sys, LINEAR_MATRIX_PARAM_NAME))
    bilinear_matrix_param = unwrap(getproperty(sys, BILINEAR_MATRIX_PARAM_NAME))
    diffcache = unwrap(getproperty(sys, DIFFCACHE_PARAM_NAME))

    floatT = calculate_float_type(op, typeof(op))
    _u0_eltype = something(u0_eltype, floatT)

    guess = copy(guesses(sys))
    guess[linear_matrix_param] = fill(NaN, size(A))
    guess[bilinear_matrix_param] = fill(NaN, size(B))
    @set! sys.guesses = guess
    defs = copy(defaults(sys))
    defs[linear_matrix_param] = A
    defs[bilinear_matrix_param] = B
    cachelen = jac ? length(x2jac) : length(x2)
    defs[diffcache] = DiffCache(zeros(DiffEqBase.value(_u0_eltype), cachelen))
    @set! sys.defaults = defs

    f, u0, p = process_SciMLProblem(SemilinearODEFunction{iip, spec}, sys, op;
        t = tspan !== nothing ? tspan[1] : tspan, expression, check_compatibility,
        semiquadratic_form, semiquadratic_jacobian, jac, sparse, u0_eltype, kwargs...)

    kwargs = process_kwargs(
        sys; expression, callback, kwargs...)

    ptype = getmetadata(sys, ProblemTypeCtx, StandardODEProblem())
    args = (; f, u0, tspan, p)
    maybe_codegen_scimlproblem(expression, SplitODEProblem{iip}, args; kwargs...)
end

function add_semilinear_parameters(sys::System)
    m = length(equations(sys))
    n = length(unknowns(sys))
    linear_matrix_param = get_linear_matrix_param((m, n))
    bilinear_matrix_param = get_bilinear_matrix_param((m, (n^2 + n) ÷ 2))
    @assert !is_parameter(sys, linear_matrix_param)
    sys = with_additional_constant_parameter(sys, linear_matrix_param)
    @assert !is_parameter(sys, bilinear_matrix_param)
    sys = with_additional_constant_parameter(sys, bilinear_matrix_param)
    @assert !is_parameter(sys, get_diffcache_param(Float64))
    diffcache = get_diffcache_param(Float64)
    sys = with_additional_nonnumeric_parameter(sys, diffcache)
    var_to_name = copy(get_var_to_name(sys))
    var_to_name[LINEAR_MATRIX_PARAM_NAME] = linear_matrix_param
    var_to_name[BILINEAR_MATRIX_PARAM_NAME] = bilinear_matrix_param
    var_to_name[DIFFCACHE_PARAM_NAME] = diffcache
    @set! sys.var_to_name = var_to_name
    if get_parent(sys) !== nothing
        @set! sys.parent = add_semilinear_parameters(get_parent(sys))
    end
    return sys
end

function check_compatible_system(
        T::Union{Type{ODEFunction}, Type{ODEProblem}, Type{DAEFunction},
            Type{DAEProblem}, Type{SteadyStateProblem}, Type{SemilinearODEFunction},
            Type{SemilinearODEProblem}},
        sys::System)
    check_time_dependent(sys, T)
    check_not_dde(sys)
    check_no_cost(sys, T)
    check_no_constraints(sys, T)
    check_no_jumps(sys, T)
    check_no_noise(sys, T)
    check_is_continuous(sys, T)
end
