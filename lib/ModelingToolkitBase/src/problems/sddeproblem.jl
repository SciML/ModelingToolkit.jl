"""$(function_docstring(SDDEFunction, true, Symbol[]))"""
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
`SciMLFunctionOptions`, bypassing the `kwargs...` wrapper above.
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

"""$(problem_docstring(SciMLBase.SDDEProblem, SDDEFunction, true))"""
@fallback_iip_specialize function SciMLBase.SDDEProblem{iip, spec}(
        sys::System, op, tspan;
        callback = nothing, check_length = true, checkbounds = false,
        eval_expression = false, eval_module = @__MODULE__, check_compatibility = true,
        u0_constructor = identity, sparse = false, sparsenoise = sparse,
        expression = Val{false}, seed = missing, constant_lags = missing, kwargs...
    ) where {iip, spec}
    fn_opts = SciMLFunctionOptions(;
        t = tspan !== nothing ? tspan[1] : tspan, eval_expression, eval_module,
        checkbounds, check_compatibility, sparse, expression, kwargs...
    )
    opts = SciMLProblemOptions(
        sys;
        fn_opts, check_length, symbolic_u0 = true, u0_constructor,
        build_initializeprob = supports_initialization(sys),
        time_dependent_init = is_time_dependent(sys),
        circular_dependency_max_cycle_length = length(all_symbols(sys)),
        kwargs...
    )
    return SDDEProblem{iip, spec}(
        sys, op, tspan, opts; callback, sparsenoise, seed, constant_lags, kwargs...
    )
end

"""
    SciMLBase.SDDEProblem{iip, spec}(sys::System, op, tspan, opts::SciMLProblemOptions; callback = nothing, sparsenoise = opts.fn_opts.sparse, kwargs...)

Public entry point that builds an `SDDEProblem` directly from a pre-assembled
[`SciMLProblemOptions`](@ref), bypassing the `kwargs...` wrapper above.

`seed`/`constant_lags` (default `missing`, meaning "not explicitly provided by the caller")
are explicit keywords here — not part of `kwargs` — since they're relevant only to the final
`SciMLBase.SDDEProblem` construction, not the inner `SDDEFunction` build; the opts-accepting
`SDDEFunction` method has no generic `kwargs...` sink to harmlessly absorb them the way its
keyword-based wrapper does.
"""
function SciMLBase.SDDEProblem{iip, spec}(
        sys::System, op, tspan, opts::SciMLProblemOptions{E};
        callback = nothing, sparsenoise = opts.fn_opts.sparse, seed = missing,
        constant_lags = missing, kwargs...
    ) where {iip, spec, E}
    check_complete(sys, SDDEProblem)
    opts.fn_opts.check_compatibility && check_compatible_system(SDDEProblem, sys)

    opts = maybe_derive_t_from_tspan(opts, tspan)

    _iip = resolve_iip(iip, op)
    f, u0,
        p = process_SciMLProblem(
        SDDEFunction{_iip, spec}, sys, op, opts; options_struct = Val(true), kwargs...
    )

    (; u0_constructor) = opts
    (; eval_expression, eval_module) = opts.fn_opts.codegen
    checkbounds = opts.fn_opts.codegen.codegen.checkbounds
    h = generate_history(
        sys, u0,
        GeneratedFunctionOptions(;
            expression = Val{E}, wrap_gfw = Val{true}, eval_expression, eval_module,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; checkbounds)
        )
    )

    if E
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

    if E
        g = :(f.g)
    else
        g = f.g
    end
    args = (; f, g, u0, h, tspan, p)
    seed_kw = seed === missing ? (;) : (; seed)
    constant_lags = resolve_constant_lags(sys, constant_lags, p)
    constant_lags_kw = constant_lags === missing ? (;) : (; constant_lags)
    kwargs = (; noise, noise_rate_prototype, seed_kw..., constant_lags_kw..., kwargs...)

    return maybe_codegen_scimlproblem(Val{E}, SDDEProblem{_iip}, args; kwargs...)
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
