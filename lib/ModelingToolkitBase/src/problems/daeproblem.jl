"""
    generate_DAENLStepData(sys, u0, p, mm, nlstep_compile, nlstep_scc; jac = false)

Generate the NLStep data for fully implicit DAE solvers. This is a stub that throws an
error if called without ModelingToolkit loaded. The actual implementation is provided by
ModelingToolkit when it is loaded.

The result is the same `SciMLBase.ODENLStepData` the mass-matrix path produces, since for
`0 = F(du, u, p, t)` both arguments of `F` are affine in the stage unknown and the six
hooks parametrize the stage system as
`g(z, p') = F(γ₁ z + outer_tmp, γ₂ z + inner_tmp, p, c)`. `γ₃` is unused and is fixed to
`1`; it is kept in the tuple only so that the setter arity matches the mass-matrix path.

When `jac = true`, the analytic Jacobian of the teared inner nonlinear system is
generated symbolically and attached to `nlprob.f.jac`, so NonlinearSolve does not
have to recompute it via AD/FD on every Newton iteration.
"""
function generate_DAENLStepData(sys, u0, p, mm, nlstep_compile, nlstep_scc; jac = false)
    error(
        """
        `nlstep=true` requires ModelingToolkit.jl to be loaded.
        Please add `using ModelingToolkit` to your code before creating a DAEProblem with `nlstep=true`.
        """
    )
end

"""$(function_docstring(DAEFunction, true, [:jac, :tgrad]))"""
@fallback_iip_specialize function SciMLBase.DAEFunction{iip, spec}(
        sys::System; u0 = nothing, p = nothing, tgrad = false, jac = false,
        t = nothing, eval_expression = false, eval_module = @__MODULE__,
        sparse = false,
        steady_state = false, checkbounds = false, sparsity = false,
        analytic = nothing,
        simplify = false, initialization_data = nothing,
        expression = Val{false}, check_compatibility = true, nlstep = false,
        nlstep_compile = true, nlstep_scc = false, optimize = nothing,
        compiler_options::CompilerOptions = CompilerOptions(), kwargs...
    ) where {iip, spec}
    opts = SciMLFunctionOptions(;
        u0, p, t, jac, tgrad, sparse, sparsity, analytic, simplify, initialization_data,
        expression, check_compatibility, eval_expression, eval_module, compiler_options,
        checkbounds, optimize, kwargs...,
    )
    return DAEFunction{iip, spec}(
        sys, opts; steady_state, nlstep, nlstep_compile, nlstep_scc
    )
end

"""
    SciMLBase.DAEFunction{iip, spec}(sys::System, opts::SciMLFunctionOptions; kwargs...)

Public entry point that builds a `DAEFunction` directly from a pre-assembled
[`SciMLFunctionOptions`](@ref), bypassing the `kwargs...` wrapper above.
"""
function SciMLBase.DAEFunction{iip, spec}(
        sys::System, opts::SciMLFunctionOptions{E};
        steady_state::Bool = false, nlstep::Bool = false, nlstep_compile::Bool = true,
        nlstep_scc::Bool = false
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

    if nlstep
        dae_nlstep = generate_DAENLStepData(
            sys, u0, p, calculate_massmatrix(sys), nlstep_compile, nlstep_scc; jac
        )
    else
        dae_nlstep = nothing
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
        nlstep_data = dae_nlstep,
    )
    args = (; f)

    return maybe_codegen_scimlfn(Val{E}, DAEFunction{iip, spec}, args; kwargs...)
end

"""$(problem_docstring(SciMLBase.DAEProblem, DAEFunction, true))"""
@fallback_iip_specialize function SciMLBase.DAEProblem{iip, spec}(
        sys::System, op, tspan;
        callback = nothing, check_length = true, eval_expression = false,
        eval_module = @__MODULE__, check_compatibility = true,
        expression = Val{false}, kwargs...
    ) where {iip, spec}
    check_complete(sys, DAEProblem)
    check_compatibility && check_compatible_system(DAEProblem, sys)

    _iip = resolve_iip(iip, op)
    f, du0,
        u0,
        p = process_SciMLProblem(
        DAEFunction{_iip, spec}, sys, op;
        t = tspan !== nothing ? tspan[1] : tspan, check_length, eval_expression,
        eval_module, check_compatibility, implicit_dae = true, expression, kwargs...
    )

    kwargs = process_kwargs(
        sys; expression, callback, eval_expression, eval_module,
        op, tspan, kwargs...
    )

    diffvars = collect_differential_variables(sys)
    sts = unknowns(sys)
    differential_vars = map(Base.Fix2(in, diffvars), sts)

    args = (; f, du0, u0, tspan, p)
    kwargs = (; differential_vars, kwargs...)

    return maybe_codegen_scimlproblem(expression, DAEProblem{_iip}, args; kwargs...)
end
