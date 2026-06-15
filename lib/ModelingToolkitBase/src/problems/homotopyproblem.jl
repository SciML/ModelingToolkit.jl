"""
    SciMLBase.HomotopyProblem(sys::System, op; λspan = (0.0, 1.0), kwargs...)

Build a [`SciMLBase.HomotopyProblem`](@ref) from a `System` whose equations
contain Modelica `homotopy(actual, simplified)` nodes (Modelica spec 3.7.4.2).

The standard nonlinear residual is built from `sys` (and used only to obtain
`u0`, `p`, the observed function, and the residual prototype), then the residual
is regenerated with every `homotopy(actual, simplified)` replaced by the convex
blend `(1 - λ)*simplified + λ*actual` and compiled as `f(u, p, λ)`. `λ` is an
explicit trailing argument — it is never added to the system's parameters, and
`p` passes through untouched — so the result solves by natural-parameter
continuation (`NonlinearSolveBase.HomotopySweep`), sweeping `λ` across `λspan`
(default `(0.0, 1.0)`, i.e. from `simplified` to `actual`).

`sys` must be `complete`d and must contain at least one `homotopy` node (use
[`SciMLBase.NonlinearProblem`](@ref) otherwise, or `AbstractNonlinearProblem` to
select automatically). The jacobian/sparsity of the standard build are
deliberately dropped: they encode the `λ = 1` (opaque-`actual`) system and would
be wrong mid-sweep.

!!! note

    `expression = Val{true}` (codegen-to-`Expr`) is not yet supported for the
    homotopy constructor; this can be added in a future PR.
"""
function SciMLBase.HomotopyProblem(
        sys::System, op;
        expression = Val{false}, λspan = (0.0, 1.0),
        check_length = true, check_compatibility = true,
        eval_expression = false, eval_module = @__MODULE__,
        checkbounds = false, cse = true, kwargs...
    )
    if expression !== Val{false}
        throw(
            ArgumentError(
                "`HomotopyProblem(sys, op)` does not yet support " *
                    "`expression = Val{true}`; build the problem directly " *
                    "(the default `expression = Val{false}`)."
            )
        )
    end
    check_complete(sys, SciMLBase.HomotopyProblem)
    if is_time_dependent(sys)
        sys = NonlinearSystem(sys)
    end
    check_compatibility && check_compatible_system(SciMLBase.NonlinearProblem, sys)
    if !has_any_homotopy(sys)
        throw(
            ArgumentError(
                "`HomotopyProblem` requires a system containing " *
                    "`homotopy(actual, simplified)` nodes; none were found. Use " *
                    "`NonlinearProblem` for systems without `homotopy`."
            )
        )
    end

    _iip = resolve_iip(Both, op)
    f, u0,
        p = process_SciMLProblem(
        SciMLBase.NonlinearFunction{_iip}, sys, op;
        check_length, check_compatibility, expression,
        eval_expression, eval_module, checkbounds, cse, kwargs...
    )

    # Swap the opaque-`actual` residual for the homotopy-swept `f(u, p, λ)`. The
    # observed function and residual prototype carry over; `initialization_data`,
    # jacobian, and sparsity are deliberately not carried (the latter two encode
    # the `λ = 1` system and would be wrong mid-sweep).
    shadow, λ = lower_homotopy(sys)
    hf = generate_homotopy_residual(
        shadow, λ; eval_expression, eval_module, checkbounds, cse
    )
    swept_f = SciMLBase.NonlinearFunction{_iip}(
        hf; sys = f.sys, observed = f.observed, resid_prototype = f.resid_prototype
    )

    kwargs = process_kwargs(sys; kwargs...)
    args = (; f = swept_f, u0, p)
    return maybe_codegen_scimlproblem(
        expression, SciMLBase.HomotopyProblem{_iip}, args; λspan, kwargs...
    )
end
