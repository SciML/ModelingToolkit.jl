@fallback_iip_specialize function SciMLBase.SDEFunction{iip, spec}(
        sys::System; u0 = nothing, p = nothing, tgrad = false, jac = false,
        t = nothing, eval_expression = false, eval_module = @__MODULE__,
        sparse = false,
        steady_state = false, checkbounds = false, sparsity = false,
        analytic = nothing,
        simplify = false, initialization_data = nothing,
        check_compatibility = true, expression = Val{false}, optimize = nothing,
        compiler_options::CompilerOptions = CompilerOptions(), kwargs...
    ) where {iip, spec}
    opts = SciMLFunctionOptions(;
        u0, p, t, jac, tgrad, sparse, sparsity, analytic, simplify, initialization_data,
        expression, check_compatibility, eval_expression, eval_module, compiler_options,
        checkbounds, optimize, kwargs...,
    )
    return SDEFunction{iip, spec}(sys, opts; steady_state)
end

"""
    SciMLBase.SDEFunction{iip, spec}(sys::System, opts::SciMLFunctionOptions; kwargs...)

Public entry point that builds an `SDEFunction` directly from a pre-assembled
[`SciMLFunctionOptions`](@ref), bypassing the `kwargs...` wrapper above.
"""
function SciMLBase.SDEFunction{iip, spec}(
        sys::System, opts::SciMLFunctionOptions{E};
        steady_state::Bool = false
    ) where {iip, spec, E}
    check_complete(sys, SDEFunction)
    opts.check_compatibility && check_compatible_system(SDEFunction, sys)

    (; u0, p, t, jac, tgrad, sparse, sparsity, analytic, simplify, initialization_data) = opts
    codegen_opts = opts.codegen

    f = generate_rhs(sys, codegen_opts)
    g = generate_diffusion_function(sys, codegen_opts)

    if spec === SciMLBase.FunctionWrapperSpecialize && iip
        if u0 === nothing || p === nothing || t === nothing
            error("u0, p, and t must be specified for FunctionWrapperSpecialize on SDEFunction.")
        end
        if E
            f = :($(SciMLBase.wrapfun_iip)($f, ($u0, $u0, $p, $t)))
        else
            f = SciMLBase.wrapfun_iip(f, (u0, u0, p, t))
        end
    end

    if tgrad
        _tgrad = generate_tgrad(sys, codegen_opts; simplify)
    else
        _tgrad = nothing
    end

    if jac
        _jac = generate_jacobian(sys, codegen_opts; simplify, sparse)
    else
        _jac = nothing
    end

    M = calculate_massmatrix(sys)
    _M = concrete_massmatrix(M; sparse, u0)

    observedfun = ObservedFunctionCache(sys, codegen_opts; steady_state)

    _W_sparsity = W_sparsity(sys)
    W_prototype = calculate_W_prototype(_W_sparsity; u0, sparse)

    kwargs = (;
        sys = sys,
        jac = _jac,
        tgrad = _tgrad,
        mass_matrix = _M,
        jac_prototype = W_prototype,
        observed = observedfun,
        sparsity = sparsity ? _W_sparsity : nothing,
        analytic = analytic,
        initialization_data,
    )
    args = (; f, g)

    return maybe_codegen_scimlfn(Val{E}, SDEFunction{iip, spec}, args; kwargs...)
end

@fallback_iip_specialize function SciMLBase.SDEProblem{iip, spec}(
        sys::System, op, tspan;
        callback = nothing, check_length = true, eval_expression = false,
        eval_module = @__MODULE__, check_compatibility = true, sparse = false,
        sparsenoise = sparse, expression = Val{false}, _skip_events = false,
        _skip_tstops = false,
        noise = missing, noise_rate_prototype = missing, seed = missing, kwargs...
    ) where {iip, spec}
    fn_opts = SciMLFunctionOptions(;
        t = tspan !== nothing ? tspan[1] : tspan, eval_expression, eval_module,
        check_compatibility, sparse, expression, kwargs...
    )
    opts = SciMLProblemOptions(
        sys;
        fn_opts, check_length, build_initializeprob = supports_initialization(sys),
        time_dependent_init = is_time_dependent(sys),
        circular_dependency_max_cycle_length = length(all_symbols(sys)),
        kwargs...
    )
    return SDEProblem{iip, spec}(
        sys, op, tspan, opts; callback, sparsenoise, _skip_events, _skip_tstops, noise,
        noise_rate_prototype, seed, kwargs...
    )
end

"""
    SciMLBase.SDEProblem{iip, spec}(sys::System, op, tspan, opts::SciMLProblemOptions; kwargs...)

Public entry point that builds an `SDEProblem` directly from a pre-assembled
[`SciMLProblemOptions`](@ref), bypassing the `kwargs...` wrapper above.

`noise`/`noise_rate_prototype`/`seed` (default `missing`, meaning "not explicitly provided
by the caller") are explicit keywords here — not part of `kwargs` — since they're relevant
only to this outer `SDEProblem` construction (`seed` reaches the final
`SciMLBase.SDEProblem` call only), not the inner `SDEFunction` build; the opts-accepting
`SDEFunction` method has no generic `kwargs...` sink to harmlessly absorb them the way its
keyword-based wrapper does. `_skip_events`/`_skip_tstops` are likewise explicit since they're
relevant only to `process_kwargs` below.
"""
function SciMLBase.SDEProblem{iip, spec}(
        sys::System, op, tspan, opts::SciMLProblemOptions{E};
        callback = nothing, sparsenoise = opts.fn_opts.sparse, _skip_events = false,
        _skip_tstops = false,
        noise = missing, noise_rate_prototype = missing, seed = missing, kwargs...
    ) where {iip, spec, E}
    check_complete(sys, SDEProblem)
    opts.fn_opts.check_compatibility && check_compatible_system(SDEProblem, sys)

    opts = maybe_derive_t_from_tspan(opts, tspan)

    _iip = resolve_iip(iip, op)
    f, u0,
        p = process_SciMLProblem(
        SDEFunction{_iip, spec}, sys, op, opts; options_struct = Val(true), kwargs...
    )

    # Only calculate noise and noise_rate_prototype if not provided by user
    if noise === missing && noise_rate_prototype === missing
        noise, noise_rate_prototype = calculate_noise_and_rate_prototype(sys, u0; sparsenoise)
    elseif noise === missing
        noise, _ = calculate_noise_and_rate_prototype(sys, u0; sparsenoise)
    elseif noise_rate_prototype === missing
        _, noise_rate_prototype = calculate_noise_and_rate_prototype(sys, u0; sparsenoise)
    end

    (; eval_expression, eval_module) = opts.fn_opts.codegen
    kwargs = process_kwargs(
        sys; expression = Val{E}, callback, eval_expression, eval_module,
        op, _skip_events, _skip_tstops, tspan, kwargs...
    )

    args = (; f, u0, tspan, p)
    seed_kw = seed === missing ? (;) : (; seed)
    kwargs = (; noise, noise_rate_prototype, seed_kw..., kwargs...)

    return maybe_codegen_scimlproblem(Val{E}, SDEProblem{_iip}, args; kwargs...)
end

function check_compatible_system(T::Union{Type{SDEFunction}, Type{SDEProblem}}, sys::System)
    check_time_dependent(sys, T)
    check_not_dde(sys)
    check_no_cost(sys, T)
    check_no_constraints(sys, T)
    check_no_jumps(sys, T)
    check_has_noise(sys, T)
    return check_is_continuous(sys, T)
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
        noise = __default_wiener_process()
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

__default_wiener_process() = __default_wiener_process(nothing)

function __default_wiener_process(_)
    error(
        """
        Generating code for this problem requires loading DiffEqNoiseProcess.jl. Please run
        `import DiffEqNoiseProcess` to proceed.
        """
    )
end
