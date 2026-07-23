function SciMLBase.IntervalNonlinearFunction(
        sys::System; u0 = nothing, p = nothing, t = nothing, eval_expression = false,
        eval_module = @__MODULE__, expression = Val{false}, checkbounds = false,
        analytic = nothing, initialization_data = nothing,
        check_compatibility = true, optimize = nothing,
        compiler_options::CompilerOptions = CompilerOptions(), kwargs...
    )
    opts = SciMLFunctionOptions(;
        u0, p, t, analytic, initialization_data,
        expression, check_compatibility, eval_expression, eval_module, compiler_options,
        checkbounds, optimize, kwargs...,
    )
    return IntervalNonlinearFunction(sys, opts)
end

"""
    SciMLBase.IntervalNonlinearFunction(sys::System, opts::SciMLFunctionOptions)

Public entry point that builds an `IntervalNonlinearFunction` directly from a pre-assembled
[`SciMLFunctionOptions`](@ref), bypassing the `kwargs...` wrapper above.
"""
function SciMLBase.IntervalNonlinearFunction(sys::System, opts::SciMLFunctionOptions{E}) where {E}
    check_complete(sys, IntervalNonlinearFunction)
    opts.check_compatibility && check_compatible_system(IntervalNonlinearFunction, sys)

    (; u0, p, analytic, initialization_data) = opts
    codegen_opts = opts.codegen

    f = generate_rhs(sys, codegen_opts; scalar = true)

    observedfun = ObservedFunctionCache(sys, codegen_opts)

    args = (; f)
    kwargs = (;
        sys = sys,
        observed = observedfun,
        analytic = analytic,
        initialization_data,
    )

    return maybe_codegen_scimlfn(
        Val{E}, IntervalNonlinearFunction{false}, args; kwargs...
    )
end

function SciMLBase.IntervalNonlinearProblem(
        sys::System, uspan::NTuple{2}, parammap = SciMLBase.NullParameters();
        check_compatibility = true, expression = Val{false}, kwargs...
    )
    check_complete(sys, IntervalNonlinearProblem)
    check_compatibility && check_compatible_system(IntervalNonlinearProblem, sys)

    u0map = unknowns(sys) .=> uspan[1]
    op = anydict([unknowns(sys)[1] => uspan[1]])
    merge!(op, to_varmap(parammap, parameters(sys)))
    f, u0,
        p = process_SciMLProblem(
        IntervalNonlinearFunction, sys, op;
        check_compatibility, expression, kwargs...
    )

    kwargs = process_kwargs(sys; kwargs...)
    args = (; f, uspan, p)
    return maybe_codegen_scimlproblem(expression, IntervalNonlinearProblem, args; kwargs...)
end

function check_compatible_system(
        T::Union{Type{IntervalNonlinearFunction}, Type{IntervalNonlinearProblem}}, sys::System
    )
    check_time_independent(sys, T)
    if !isone(length(unknowns(sys)))
        throw(
            SystemCompatibilityError(
                """
                    `$T` requires a system with a single unknown. Found `$(unknowns(sys))`.
                """
            )
        )
    end
    if !isone(length(equations(sys)))
        throw(
            SystemCompatibilityError(
                """
                    `$T` requires a system with a single equation. Found `$(equations(sys))`.
                """
            )
        )
    end
    check_not_dde(sys)
    check_no_cost(sys, T)
    check_no_constraints(sys, T)
    check_no_jumps(sys, T)
    return check_no_noise(sys, T)
end
