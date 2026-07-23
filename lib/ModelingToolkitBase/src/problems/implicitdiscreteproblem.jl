@fallback_iip_specialize function SciMLBase.ImplicitDiscreteFunction{iip, spec}(
        sys::System; u0 = nothing, p = nothing, t = nothing, eval_expression = false,
        eval_module = @__MODULE__, expression = Val{false},
        checkbounds = false, analytic = nothing, simplify = false,
        initialization_data = nothing, check_compatibility = true,
        optimize = nothing, compiler_options::CompilerOptions = CompilerOptions(), kwargs...
    ) where {
        iip, spec,
    }
    opts = SciMLFunctionOptions(;
        u0, p, t, analytic, simplify, initialization_data,
        expression, check_compatibility, eval_expression, eval_module, compiler_options,
        checkbounds, optimize, kwargs...,
    )
    return ImplicitDiscreteFunction{iip, spec}(sys, opts)
end

"""
    SciMLBase.ImplicitDiscreteFunction{iip, spec}(sys::System, opts::SciMLFunctionOptions)

Public entry point that builds an `ImplicitDiscreteFunction` directly from a pre-assembled
[`SciMLFunctionOptions`](@ref), bypassing the `kwargs...` wrapper above.
"""
function SciMLBase.ImplicitDiscreteFunction{iip, spec}(
        sys::System, opts::SciMLFunctionOptions{E}
    ) where {iip, spec, E}
    check_complete(sys, ImplicitDiscreteFunction)
    opts.check_compatibility && check_compatible_system(ImplicitDiscreteFunction, sys)

    iv = get_iv(sys)
    dvs = unknowns(sys)
    (; u0, p, t, analytic, initialization_data) = opts
    codegen_opts = opts.codegen

    f = generate_rhs(sys, codegen_opts; implicit_dae = true, override_discrete = true)

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

    observedfun = ObservedFunctionCache(sys, codegen_opts)

    args = (; f)
    kwargs = (;
        sys = sys,
        observed = observedfun,
        analytic = analytic,
        initialization_data,
        resid_prototype,
    )

    return maybe_codegen_scimlfn(
        Val{E}, ImplicitDiscreteFunction{iip, spec}, args; kwargs...
    )
end

@fallback_iip_specialize function SciMLBase.ImplicitDiscreteProblem{iip, spec}(
        sys::System, op, tspan;
        check_compatibility = true, expression = Val{false}, kwargs...
    ) where {iip, spec}
    check_complete(sys, ImplicitDiscreteProblem)
    check_compatibility && check_compatible_system(ImplicitDiscreteProblem, sys)

    _iip = resolve_iip(iip, op)
    dvs = unknowns(sys)
    op = to_varmap(op, dvs)
    add_toterms!(op; replace = true)
    f, u0,
        p = process_SciMLProblem(
        ImplicitDiscreteFunction{_iip, spec}, sys, op;
        t = tspan !== nothing ? tspan[1] : tspan, check_compatibility,
        expression, kwargs...
    )

    kwargs = process_kwargs(sys; kwargs...)
    args = (; f, u0, tspan, p)
    return maybe_codegen_scimlproblem(
        expression, ImplicitDiscreteProblem{_iip}, args; kwargs...
    )
end

function check_compatible_system(
        T::Union{Type{ImplicitDiscreteFunction}, Type{ImplicitDiscreteProblem}}, sys::System
    )
    check_time_dependent(sys, T)
    check_not_dde(sys)
    check_no_cost(sys, T)
    check_no_constraints(sys, T)
    check_no_jumps(sys, T)
    check_no_noise(sys, T)
    return check_is_discrete(sys, T)
end
