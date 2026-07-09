# Backwards-compatibility keyword-argument wrappers for the code-generation entry points
# defined in `ModelingToolkit` (as opposed to `ModelingToolkitBase`). See
# `ModelingToolkitBase`'s `systems/codegen_compat.jl` for the rationale: the primary methods
# take a positional `GeneratedFunctionOptions` plus strictly-typed function-specific keyword
# arguments; these wrappers accept the historical loose keyword arguments and forward.

function generate_semiquadratic_functions(
        sys::System, A, B, C; stiff_linear = true, stiff_quadratic = false,
        stiff_nonlinear = false, expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__, kwargs...
    )
    return generate_semiquadratic_functions(
        sys, A, B, C,
        GeneratedFunctionOptions(;
            expression, wrap_gfw, eval_expression, eval_module,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; kwargs...)
        );
        stiff_linear, stiff_quadratic, stiff_nonlinear
    )
end

function generate_semiquadratic_jacobian(
        sys::System, A, B, C, Cjac; sparse = false, stiff_linear = true,
        stiff_quadratic = false, stiff_nonlinear = false, expression = Val{true},
        wrap_gfw = Val{false}, eval_expression = false, eval_module = @__MODULE__, kwargs...
    )
    return generate_semiquadratic_jacobian(
        sys, A, B, C, Cjac,
        GeneratedFunctionOptions(;
            expression, wrap_gfw, eval_expression, eval_module,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; kwargs...)
        );
        sparse, stiff_linear, stiff_quadratic, stiff_nonlinear
    )
end
