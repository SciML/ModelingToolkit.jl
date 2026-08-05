@fallback_iip_specialize function SciMLBase.DAEFunction{iip, spec}(
        sys::System; u0 = nothing, p = nothing, tgrad = false, jac = false,
        t = nothing, eval_expression = false, eval_module = @__MODULE__,
        sparse = false,
        steady_state = false, checkbounds = false, sparsity = false,
        analytic = nothing,
        simplify = false, initialization_data = nothing,
        expression = Val{false}, check_compatibility = true, optimize = nothing,
        compiler_options::CompilerOptions = CompilerOptions(), kwargs...
    ) where {iip, spec}
    opts = SciMLFunctionOptions(;
        u0, p, t, jac, tgrad, sparse, sparsity, analytic, simplify, initialization_data,
        expression, check_compatibility, eval_expression, eval_module, compiler_options,
        checkbounds, optimize, kwargs...,
    )
    return DAEFunction{iip, spec}(sys, opts; steady_state)
end

"""
    SciMLBase.DAEFunction{iip, spec}(sys::System, opts::SciMLFunctionOptions; kwargs...)

Public entry point that builds a `DAEFunction` directly from a pre-assembled
[`SciMLFunctionOptions`](@ref), bypassing the `kwargs...` wrapper above.
"""
function SciMLBase.DAEFunction{iip, spec}(
        sys::System, opts::SciMLFunctionOptions{E};
        steady_state::Bool = false
    ) where {iip, spec, E}
    check_complete(sys, DAEFunction)
    opts.check_compatibility && check_compatible_system(DAEFunction, sys)

    (; u0, p, t, jac, tgrad, sparse, analytic, simplify, initialization_data) = opts
    codegen_opts = opts.codegen

    f = generate_rhs(sys, codegen_opts; implicit_dae = true)

    if spec === SciMLBase.FunctionWrapperSpecialize && iip
        if u0 === nothing || p === nothing || t === nothing
            error("u0, p, and t must be specified for FunctionWrapperSpecialize on ODEFunction.")
        end
        if E
            f = :($(SciMLBase.wrapfun_iip)($f, ($u0, $u0, $u0, $p, $t)))
        else
            f = SciMLBase.wrapfun_iip(f, (u0, u0, u0, p, t))
        end
    end

    if jac
        _jac = generate_dae_jacobian(sys, codegen_opts; simplify, sparse)
    else
        _jac = nothing
    end

    observedfun = ObservedFunctionCache(sys, codegen_opts; steady_state)

    jac_prototype = if sparse
        uElType = u0 === nothing ? Float64 : eltype(u0)
        if jac
            J1 = calculate_jacobian(sys, sparse = sparse)
            derivatives = Differential(get_iv(sys)).(unknowns(sys))
            J2 = calculate_jacobian(sys; sparse = sparse, dvs = derivatives)
            similar(J1 + J2, uElType)
        else
            similar(jacobian_dae_sparsity(sys), uElType)
        end
    else
        nothing
    end

    kwargs = (;
        sys = sys,
        jac = _jac,
        jac_prototype = jac_prototype,
        observed = observedfun,
        analytic = analytic,
        initialization_data,
    )
    args = (; f)

    return maybe_codegen_scimlfn(Val{E}, DAEFunction{iip, spec}, args; kwargs...)
end

@fallback_iip_specialize function SciMLBase.DAEProblem{iip, spec}(
        sys::System, op, tspan;
        callback = nothing, check_length = true, eval_expression = false,
        eval_module = @__MODULE__, check_compatibility = true,
        expression = Val{false}, kwargs...
    ) where {iip, spec}
    fn_opts = SciMLFunctionOptions(;
        t = tspan !== nothing ? tspan[1] : tspan, eval_expression, eval_module,
        check_compatibility, expression, kwargs...
    )
    opts = SciMLProblemOptions(
        sys;
        fn_opts, check_length, implicit_dae = true,
        build_initializeprob = supports_initialization(sys),
        time_dependent_init = is_time_dependent(sys),
        circular_dependency_max_cycle_length = length(all_symbols(sys)),
        kwargs...
    )
    return DAEProblem{iip, spec}(sys, op, tspan, opts; callback, kwargs...)
end

"""
    SciMLBase.DAEProblem{iip, spec}(sys::System, op, tspan, opts::SciMLProblemOptions; callback = nothing, kwargs...)

Public entry point that builds a `DAEProblem` directly from a pre-assembled
[`SciMLProblemOptions`](@ref), bypassing the `kwargs...` wrapper above.
"""
function SciMLBase.DAEProblem{iip, spec}(
        sys::System, op, tspan, opts::SciMLProblemOptions{E};
        callback = nothing, kwargs...
    ) where {iip, spec, E}
    check_complete(sys, DAEProblem)
    opts.fn_opts.check_compatibility && check_compatible_system(DAEProblem, sys)

    opts = maybe_derive_t_from_tspan(opts, tspan)

    _iip = resolve_iip(iip, op)
    f, du0,
        u0,
        p = process_SciMLProblem(
        DAEFunction{_iip, spec}, sys, op, opts; options_struct = Val(true), kwargs...
    )

    (; eval_expression, eval_module) = opts.fn_opts.codegen
    kwargs = process_kwargs(
        sys; expression = Val{E}, callback, eval_expression, eval_module,
        op, tspan, kwargs...
    )

    diffvars = collect_differential_variables(sys)
    sts = unknowns(sys)
    differential_vars = map(Base.Fix2(in, diffvars), sts)

    args = (; f, du0, u0, tspan, p)
    kwargs = (; differential_vars, kwargs...)

    return maybe_codegen_scimlproblem(Val{E}, DAEProblem{_iip}, args; kwargs...)
end
