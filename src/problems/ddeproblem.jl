@fallback_iip_specialize function SciMLBase.DDEFunction{iip, spec}(
        sys::System; u0 = nothing, p = nothing, eval_expression = false,
        eval_module = @__MODULE__, expression = Val{false}, checkbounds = false,
        initialization_data = nothing, cse = true, check_compatibility = true,
        sparse = false, simplify = false, analytic = nothing, kwargs...) where {iip, spec}
    check_complete(sys, DDEFunction)
    check_compatibility && check_compatible_system(DDEFunction, sys)

    f = generate_rhs(sys; expression, wrap_gfw = Val{true},
        eval_expression, eval_module, checkbounds = checkbounds, cse,
        kwargs...)

    if spec === SciMLBase.FunctionWrapperSpecialize && iip
        if u0 === nothing || p === nothing || t === nothing
            error("u0, p, and t must be specified for FunctionWrapperSpecialize on DDEFunction.")
        end
        if expression == Val{true}
            f = :($(SciMLBase.wrapfun_iip)($f, ($u0, $u0, $p, $t)))
        else
            f = SciMLBase.wrapfun_iip(f, (u0, u0, p, t))
        end
    end

    M = calculate_massmatrix(sys)
    _M = concrete_massmatrix(M; sparse, u0)

    observedfun = ObservedFunctionCache(
        sys; expression, eval_expression, eval_module, checkbounds, cse)

    kwargs = (;
        sys = sys,
        mass_matrix = _M,
        observed = observedfun,
        analytic = analytic,
        initialization_data)
    args = (; f)

    return maybe_codegen_scimlfn(expression, DDEFunction{iip, spec}, args; kwargs...)
end

@fallback_iip_specialize function SciMLBase.DDEProblem{iip, spec}(
        sys::System, op, tspan;
        callback = nothing, check_length = true, cse = true, checkbounds = false,
        eval_expression = false, eval_module = @__MODULE__, check_compatibility = true,
        u0_constructor = identity, expression = Val{false}, kwargs...) where {iip, spec}
    check_complete(sys, DDEProblem)
    check_compatibility && check_compatible_system(DDEProblem, sys)

    f, u0,
    p = process_SciMLProblem(DDEFunction{iip, spec}, sys, op;
        t = tspan !== nothing ? tspan[1] : tspan, check_length, cse, checkbounds,
        eval_expression, eval_module, check_compatibility, symbolic_u0 = true,
        expression, u0_constructor, kwargs...)

    h = generate_history(
        sys, u0; expression, wrap_gfw = Val{true}, cse, eval_expression, eval_module,
        checkbounds)

    if expression == Val{true}
        if u0 !== nothing
            u0 = :($u0_constructor($map($float, h(p, tspan[1]))))
        end
    else
        if u0 !== nothing
            u0 = u0_constructor(float.(h(p, tspan[1])))
        end
    end

    kwargs = process_kwargs(
        sys; expression, callback, eval_expression, eval_module, op, kwargs...)
    args = (; f, u0, h, tspan, p)

    return maybe_codegen_scimlproblem(expression, DDEProblem{iip}, args; kwargs...)
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
