# Backwards-compatibility keyword-argument wrappers for the code-generation entry points.
#
# The primary methods (in `codegen.jl`) take a positional `GeneratedFunctionOptions` for the
# output/compile options plus a strictly-typed set of function-specific keyword arguments.
# Each wrapper here accepts the historical loose keyword arguments, bundles the
# code-generation options into a `Symbolics.CodegenFunctionOptions` and the output/compile
# options into a `GeneratedFunctionOptions`, and forwards to the corresponding primary
# method. New code should construct a `GeneratedFunctionOptions` and call the primary method
# directly.

function generate_rhs(
        sys::System; expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__,
        compiler_options::CompilerOptions = CompilerOptions(),
        implicit_dae = false, scalar = false, override_discrete = false,
        cachesyms = nothing, extra_args::Tuple = (), kwargs...
    )
    return generate_rhs(
        sys,
        GeneratedFunctionOptions(;
            expression, wrap_gfw, eval_expression, eval_module,
            compiler_options,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; kwargs...)
        );
        implicit_dae, scalar, override_discrete, cachesyms, extra_args
    )
end

function generate_diffusion_function(
        sys::System; expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__, kwargs...
    )
    return generate_diffusion_function(
        sys,
        GeneratedFunctionOptions(;
            expression, wrap_gfw, eval_expression, eval_module,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; kwargs...)
        )
    )
end

function generate_jacobian(
        sys::System; expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__,
        compiler_options::CompilerOptions = CompilerOptions(),
        simplify = false, sparse = false, kwargs...
    )
    return generate_jacobian(
        sys,
        GeneratedFunctionOptions(;
            expression, wrap_gfw, eval_expression, eval_module,
            compiler_options,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; kwargs...)
        );
        simplify, sparse
    )
end

function generate_tgrad(
        sys::System; expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__,
        compiler_options::CompilerOptions = CompilerOptions(),
        simplify = false, kwargs...
    )
    return generate_tgrad(
        sys,
        GeneratedFunctionOptions(;
            expression, wrap_gfw, eval_expression, eval_module,
            compiler_options,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; kwargs...)
        );
        simplify
    )
end

function generate_W(
        sys::System; expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__,
        simplify = false, sparse = false, kwargs...
    )
    return generate_W(
        sys,
        GeneratedFunctionOptions(;
            expression, wrap_gfw, eval_expression, eval_module,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; kwargs...)
        );
        simplify, sparse
    )
end

function generate_dae_jacobian(
        sys::System; expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__,
        compiler_options::CompilerOptions = CompilerOptions(),
        simplify = false, sparse = false, kwargs...
    )
    return generate_dae_jacobian(
        sys,
        GeneratedFunctionOptions(;
            expression, wrap_gfw, eval_expression, eval_module,
            compiler_options,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; kwargs...)
        );
        simplify, sparse
    )
end

function generate_cost(
        sys::System; expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__, kwargs...
    )
    return generate_cost(
        sys,
        GeneratedFunctionOptions(;
            expression, wrap_gfw, eval_expression, eval_module,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; kwargs...)
        )
    )
end

function generate_cons(
        sys::System; expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__, kwargs...
    )
    return generate_cons(
        sys,
        GeneratedFunctionOptions(;
            expression, wrap_gfw, eval_expression, eval_module,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; kwargs...)
        )
    )
end

function generate_bvp_cost(
        sys::System; expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__, kwargs...
    )
    return generate_bvp_cost(
        sys,
        GeneratedFunctionOptions(;
            expression, wrap_gfw, eval_expression, eval_module,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; kwargs...)
        )
    )
end

function generate_cost_gradient(
        sys::System; expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__, simplify = false, kwargs...
    )
    return generate_cost_gradient(
        sys,
        GeneratedFunctionOptions(;
            expression, wrap_gfw, eval_expression, eval_module,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; kwargs...)
        );
        simplify
    )
end

function generate_cost_hessian(
        sys::System; expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__, simplify = false,
        sparse = false, return_sparsity = false, kwargs...
    )
    return generate_cost_hessian(
        sys,
        GeneratedFunctionOptions(;
            expression, wrap_gfw, eval_expression, eval_module,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; kwargs...)
        );
        simplify, sparse, return_sparsity
    )
end

function generate_constraint_jacobian(
        sys::System; expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__, return_sparsity = false,
        simplify = false, sparse = false, kwargs...
    )
    return generate_constraint_jacobian(
        sys,
        GeneratedFunctionOptions(;
            expression, wrap_gfw, eval_expression, eval_module,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; kwargs...)
        );
        return_sparsity, simplify, sparse
    )
end

function generate_constraint_hessian(
        sys::System; expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__, return_sparsity = false,
        simplify = false, sparse = false, kwargs...
    )
    return generate_constraint_hessian(
        sys,
        GeneratedFunctionOptions(;
            expression, wrap_gfw, eval_expression, eval_module,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; kwargs...)
        );
        return_sparsity, simplify, sparse
    )
end

function generate_control_jacobian(
        sys::AbstractSystem; expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__, simplify = false,
        sparse = false, kwargs...
    )
    return generate_control_jacobian(
        sys,
        GeneratedFunctionOptions(;
            expression, wrap_gfw, eval_expression, eval_module,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; kwargs...)
        );
        simplify, sparse
    )
end

function generate_history(
        sys::System, u0; expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__, kwargs...
    )
    return generate_history(
        sys, u0,
        GeneratedFunctionOptions(;
            expression, wrap_gfw, eval_expression, eval_module,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; kwargs...)
        )
    )
end

function generate_boundary_conditions(
        sys::System, u0, u0_idxs, t0; expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__, kwargs...
    )
    return generate_boundary_conditions(
        sys, u0, u0_idxs, t0,
        GeneratedFunctionOptions(;
            expression, wrap_gfw, eval_expression, eval_module,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; kwargs...)
        )
    )
end

for AType in (
        :AbstractMatrix, :(Diagonal{SymbolicT, Vector{SymbolicT}}),
        :(BandedMatrix{SymbolicT, Matrix{SymbolicT}}),
    )
    @eval function generate_update_A(
            sys::System, A::$AType; expression = Val{true}, wrap_gfw = Val{false},
            eval_expression = false, eval_module = @__MODULE__,
            compiler_options::CompilerOptions = CompilerOptions(), cachesyms = (), kwargs...
        )
        return generate_update_A(
            sys, A,
            GeneratedFunctionOptions(;
                expression, wrap_gfw, eval_expression, eval_module,
                compiler_options,
                codegen_function_options = Symbolics.CodegenFunctionOptions(; kwargs...)
            );
            cachesyms
        )
    end
end

function generate_update_b(
        sys::System, b::AbstractVector; expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__,
        compiler_options::CompilerOptions = CompilerOptions(), cachesyms = (), kwargs...
    )
    return generate_update_b(
        sys, b,
        GeneratedFunctionOptions(;
            expression, wrap_gfw, eval_expression, eval_module,
            compiler_options,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; kwargs...)
        );
        cachesyms
    )
end

function generate_custom_function(
        sys::AbstractSystem, exprs, dvs = unknowns(sys),
        ps = parameters(sys; initial_parameters = true);
        expression = Val{true}, eval_expression = false, eval_module = @__MODULE__,
        cachesyms::Tuple = (), kwargs...
    )
    return generate_custom_function(
        sys, exprs, dvs, ps,
        GeneratedFunctionOptions(;
            expression, eval_expression, eval_module,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; kwargs...)
        );
        cachesyms
    )
end

Base.@nospecializeinfer function build_explicit_observed_function(
        sys, @nospecialize(ts);
        inputs = nothing, disturbance_inputs = nothing, known_disturbance_inputs = nothing,
        disturbance_argument = false, expression = false, eval_expression = false,
        eval_module = @__MODULE__, output_type = Array, checkbounds = false,
        ps = parameters(sys; initial_parameters = true), return_inplace = Val(false),
        param_only = false, throw = true, wrap_delays = is_dde(sys) && !param_only,
        force_time_independent = false, compiler_options::CompilerOptions = CompilerOptions(),
        kwargs...
    )
    return build_explicit_observed_function(
        sys, ts,
        GeneratedFunctionOptions(;
            expression, eval_expression, eval_module, compiler_options,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; checkbounds, kwargs...)
        );
        inputs, disturbance_inputs, known_disturbance_inputs, disturbance_argument,
        output_type, ps, return_inplace, param_only, throw, wrap_delays, force_time_independent
    )
end

Base.@nospecializeinfer function compile_condition(
        @nospecialize(cbs::Union{AbstractCallback, Vector{<:AbstractCallback}}),
        sys, @nospecialize(dvs), @nospecialize(ps);
        eval_expression = false, eval_module = @__MODULE__, kwargs...
    )
    return compile_condition(
        cbs, sys, dvs, ps,
        GeneratedFunctionOptions(;
            eval_expression, eval_module,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; kwargs...)
        )
    )
end

Base.@nospecializeinfer function compile_explicit_affect(
        @nospecialize(aff::AffectSystem), sys;
        reset_jumps = false, eval_expression = false, eval_module = @__MODULE__, kwargs...
    )
    return compile_explicit_affect(
        aff, sys,
        GeneratedFunctionOptions(;
            eval_expression, eval_module,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; kwargs...)
        );
        reset_jumps
    )
end

function generate_control_function(
        sys::AbstractSystem, inputs = default_codegen_inputs(sys),
        disturbance_inputs = disturbances(sys);
        known_disturbance_inputs = nothing, disturbance_argument = false,
        implicit_dae = false, simplify = false, eval_expression = false,
        eval_module = @__MODULE__, split = true, kwargs...
    )
    return generate_control_function(
        sys, inputs, disturbance_inputs,
        GeneratedFunctionOptions(;
            eval_expression, eval_module,
            codegen_function_options = Symbolics.CodegenFunctionOptions(; kwargs...)
        );
        known_disturbance_inputs, disturbance_argument, implicit_dae, simplify, split
    )
end
