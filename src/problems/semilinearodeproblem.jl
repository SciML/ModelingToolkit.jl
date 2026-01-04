@fallback_iip_specialize function SemilinearODEFunction{iip, specialize}(
        sys::System; u0 = nothing, p = nothing, t = nothing,
        semiquadratic_form = nothing,
        stiff_linear = true, stiff_quadratic = false, stiff_nonlinear = false,
        eval_expression = false, eval_module = @__MODULE__,
        expression = Val{false}, sparse = false, check_compatibility = true,
        jac = false, checkbounds = false, cse = true, initialization_data = nothing,
        analytic = nothing, kwargs...
    ) where {iip, specialize}
    check_complete(sys, SemilinearODEFunction)
    check_compatibility && check_compatible_system(SemilinearODEFunction, sys)

    if semiquadratic_form === nothing
        semiquadratic_form = calculate_semiquadratic_form(sys; sparse)
        sys = add_semiquadratic_parameters(sys, semiquadratic_form...)
    end

    A, B, C = semiquadratic_form
    M = calculate_massmatrix(sys)
    _M = concrete_massmatrix(M; sparse, u0)
    dvs = unknowns(sys)

    f1,
        f2 = generate_semiquadratic_functions(
        sys, A, B, C; stiff_linear, stiff_quadratic,
        stiff_nonlinear, expression, wrap_gfw = Val{true},
        eval_expression, eval_module, kwargs...
    )

    if jac
        Cjac = (C === nothing || !stiff_nonlinear) ? nothing : Symbolics.jacobian(C, dvs)
        _jac = generate_semiquadratic_jacobian(
            sys, A, B, C, Cjac; sparse, expression,
            wrap_gfw = Val{true}, eval_expression, eval_module, kwargs...
        )
        _W_sparsity = get_semiquadratic_W_sparsity(
            sys, A, B, C, Cjac; stiff_linear, stiff_quadratic, stiff_nonlinear, mm = M
        )
        W_prototype = calculate_W_prototype(_W_sparsity; u0, sparse)
    else
        _jac = nothing
        W_prototype = nothing
    end

    observedfun = ObservedFunctionCache(
        sys; expression, steady_state = false, eval_expression, eval_module, checkbounds, cse
    )

    args = (; f1)
    kwargs = (; jac = _jac, jac_prototype = W_prototype)
    f1 = maybe_codegen_scimlfn(expression, ODEFunction{iip, specialize}, args; kwargs...)

    args = (; f1, f2)
    kwargs = (;
        sys = sys,
        jac = _jac,
        mass_matrix = _M,
        jac_prototype = W_prototype,
        observed = observedfun,
        analytic,
        initialization_data,
    )

    return maybe_codegen_scimlfn(
        expression, SplitFunction{iip, specialize}, args; kwargs...
    )
end

@fallback_iip_specialize function SemilinearODEProblem{iip, spec}(
        sys::System, op, tspan; check_compatibility = true, u0_eltype = nothing,
        expression = Val{false}, callback = nothing, sparse = false,
        stiff_linear = true, stiff_quadratic = false, stiff_nonlinear = false,
        jac = false, kwargs...
    ) where {
        iip, spec,
    }
    check_complete(sys, SemilinearODEProblem)
    check_compatibility && check_compatible_system(SemilinearODEProblem, sys)
    if sparse
        error(
            """
            The sparse form for `SemilinearODEProblem` is currently unavailable.
            """
        )
    end

    A, B, C = semiquadratic_form = calculate_semiquadratic_form(sys; sparse)
    eqs = equations(sys)
    dvs = unknowns(sys)

    sys = add_semiquadratic_parameters(sys, A, B, C)
    if A !== nothing
        linear_matrix_param = unwrap(getproperty(sys, LINEAR_MATRIX_PARAM_NAME))
    else
        linear_matrix_param = nothing
    end
    if B !== nothing
        quadratic_forms = [
            unwrap(getproperty(sys, get_quadratic_form_name(i)))
                for i in 1:length(eqs)
        ]
        diffcache_par = unwrap(getproperty(sys, DIFFCACHE_PARAM_NAME))
    else
        quadratic_forms = diffcache_par = nothing
    end

    op = to_varmap(op, dvs)
    floatT = calculate_float_type(op, typeof(op))
    _u0_eltype = something(u0_eltype, floatT)

    guess = copy(guesses(sys))
    defs = copy(initial_conditions(sys))
    if A !== nothing
        guess[linear_matrix_param] = fill(NaN, size(A))
        defs[linear_matrix_param] = A
    end
    if B !== nothing
        for (par, mat) in zip(quadratic_forms, B)
            guess[par] = fill(NaN, size(mat))
            defs[par] = mat
        end
        cachelen = jac ? length(dvs) * length(eqs) : length(dvs)
        defs[diffcache_par] = DiffCache(zeros(DiffEqBase.value(_u0_eltype), cachelen))
    end
    @set! sys.guesses = guess
    @set! sys.initial_conditions = defs

    f, u0,
        p = process_SciMLProblem(
        SemilinearODEFunction{iip, spec}, sys, op;
        t = tspan !== nothing ? tspan[1] : tspan, expression, check_compatibility,
        semiquadratic_form, sparse, u0_eltype, stiff_linear, stiff_quadratic, stiff_nonlinear, jac, kwargs...
    )

    kwargs = process_kwargs(sys; expression, callback, kwargs...)

    args = (; f, u0, tspan, p)
    maybe_codegen_scimlproblem(expression, SplitODEProblem{iip}, args; kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Add the necessary parameters for [`SemilinearODEProblem`](@ref) given the matrices
`A`, `B`, `C` returned from [`calculate_semiquadratic_form`](@ref).
"""
function add_semiquadratic_parameters(sys::System, A, B, C)
    eqs = equations(sys)
    n = length(eqs)
    var_to_name = copy(get_var_to_name(sys))
    if B !== nothing
        for i in eachindex(B)
            B[i] === nothing && continue
            par = get_quadratic_form_param((n, n), i)
            var_to_name[get_quadratic_form_name(i)] = par
            sys = with_additional_constant_parameter(sys, par)
        end
        par = get_diffcache_param(Float64)
        var_to_name[DIFFCACHE_PARAM_NAME] = par
        sys = with_additional_nonnumeric_parameter(sys, par)
    end
    if A !== nothing
        par = get_linear_matrix_param((n, n))
        var_to_name[LINEAR_MATRIX_PARAM_NAME] = par
        sys = with_additional_constant_parameter(sys, par)
    end
    @set! sys.var_to_name = var_to_name
    if get_parent(sys) !== nothing
        @set! sys.parent = add_semiquadratic_parameters(get_parent(sys), A, B, C)
    end
    return sys
end

function MTKBase.check_compatible_system(
        T::Union{Type{SemilinearODEFunction}, Type{SemilinearODEProblem}},
        sys::System
    )
    check_time_dependent(sys, T)
    check_not_dde(sys)
    check_no_cost(sys, T)
    check_no_constraints(sys, T)
    check_no_jumps(sys, T)
    check_no_noise(sys, T)
    return check_is_continuous(sys, T)
end
