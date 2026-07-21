@fallback_iip_specialize function SciMLBase.SDDEFunction{iip, spec}(
        sys::System; u0 = nothing, p = nothing, t = nothing, expression = Val{false},
        eval_expression = false, eval_module = @__MODULE__,
        checkbounds = false,
        initialization_data = nothing, check_compatibility = true,
        sparse = false, simplify = false, analytic = nothing,
        optimize = nothing, compiler_options::CompilerOptions = CompilerOptions(), kwargs...
    ) where {iip, spec}
    opts = SciMLFunctionOptions(;
        u0, p, t, sparse, analytic, simplify, initialization_data,
        expression, check_compatibility, eval_expression, eval_module, compiler_options,
        checkbounds, optimize, kwargs...,
    )
    return SDDEFunction{iip, spec}(sys, opts)
end

"""
    SciMLBase.SDDEFunction{iip, spec}(sys::System, opts::SciMLFunctionOptions)

Public entry point that builds an `SDDEFunction` directly from a pre-assembled
[`SciMLFunctionOptions`](@ref), bypassing the `kwargs...` wrapper above.
"""
function SciMLBase.SDDEFunction{iip, spec}(
        sys::System, opts::SciMLFunctionOptions{E}
    ) where {iip, spec, E}
    check_complete(sys, SDDEFunction)
    opts.check_compatibility && check_compatible_system(SDDEFunction, sys)

    (; u0, p, t, sparse, analytic, initialization_data) = opts
    codegen_opts = opts.codegen

    f = generate_rhs(sys, codegen_opts)
    g = generate_diffusion_function(sys, codegen_opts)

    if spec === SciMLBase.FunctionWrapperSpecialize && iip
        if u0 === nothing || p === nothing || t === nothing
            error("u0, p, and t must be specified for FunctionWrapperSpecialize on SDDEFunction.")
        end
        if E
            f = :($(SciMLBase.wrapfun_iip)($f, ($u0, $u0, $p, $t)))
        else
            f = SciMLBase.wrapfun_iip(f, (u0, u0, p, t))
        end
    end

    M = calculate_massmatrix(sys)
    _M = concrete_massmatrix(M; sparse, u0)

    observedfun = ObservedFunctionCache(sys, codegen_opts)

    kwargs = (;
        sys = sys,
        mass_matrix = _M,
        observed = observedfun,
        analytic = analytic,
        initialization_data,
    )
    args = (; f, g)

    return maybe_codegen_scimlfn(Val{E}, SDDEFunction{iip, spec}, args; kwargs...)
end

@fallback_iip_specialize function SciMLBase.SDDEProblem{iip, spec}(
        sys::System, op, tspan;
        callback = nothing, check_length = true, checkbounds = false,
        eval_expression = false, eval_module = @__MODULE__, check_compatibility = true,
        u0_constructor = identity, sparse = false, sparsenoise = sparse,
        expression = Val{false}, kwargs...
    ) where {iip, spec}
    check_complete(sys, SDDEProblem)
    check_compatibility && check_compatible_system(SDDEProblem, sys)

    _iip = resolve_iip(iip, op)
    f, u0,
        p = process_SciMLProblem(
        SDDEFunction{_iip, spec}, sys, op;
        t = tspan !== nothing ? tspan[1] : tspan, check_length, checkbounds,
        eval_expression, eval_module, check_compatibility, sparse, symbolic_u0 = true,
        expression, u0_constructor, kwargs...
    )

    h = generate_history(
        sys, u0,
        GeneratedFunctionOptions(;
            expression, wrap_gfw = Val{true}, eval_expression, eval_module,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; checkbounds)
        )
    )

    if expression == Val{true}
        if u0 !== nothing
            u0 = :($u0_constructor($map($float, h(p, tspan[1]))))
        end
    else
        if u0 !== nothing
            u0 = u0_constructor(float.(h(p, tspan[1])))
        end
    end

    noise, noise_rate_prototype = calculate_noise_and_rate_prototype(sys, u0; sparsenoise)
    kwargs = process_kwargs(sys; callback, eval_expression, eval_module, op, tspan, kwargs...)

    if expression == Val{true}
        g = :(f.g)
    else
        g = f.g
    end
    args = (; f, g, u0, h, tspan, p)
    kwargs = (; noise, noise_rate_prototype, kwargs...)

    return maybe_codegen_scimlproblem(expression, SDDEProblem{_iip}, args; kwargs...)
end

function check_compatible_system(
        T::Union{Type{SDDEFunction}, Type{SDDEProblem}}, sys::System
    )
    check_time_dependent(sys, T)
    check_is_dde(sys)
    check_no_cost(sys, T)
    check_no_constraints(sys, T)
    check_no_jumps(sys, T)
    check_has_noise(sys, T)
    return check_is_continuous(sys, T)
end
