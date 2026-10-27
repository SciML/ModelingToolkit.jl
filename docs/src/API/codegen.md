# Code generation utilities

These are lower-level functions that ModelingToolkit leverages to generate code for
building numerical problems.

```@docs
ModelingToolkit.generate_rhs
ModelingToolkit.generate_diffusion_function
ModelingToolkit.generate_jacobian
ModelingToolkit.generate_tgrad
ModelingToolkit.generate_W
ModelingToolkit.generate_dae_jacobian
ModelingToolkit.generate_history
ModelingToolkit.generate_boundary_conditions
ModelingToolkitBase.generate_trajectory
ModelingToolkit.generate_cost
ModelingToolkit.generate_cost_gradient
ModelingToolkit.generate_cost_hessian
ModelingToolkit.generate_cons
ModelingToolkit.generate_constraint_jacobian
ModelingToolkit.generate_constraint_hessian
ModelingToolkit.generate_multiobjective_cost
ModelingToolkit.generate_multiobjective_jacobian
ModelingToolkit.generate_multiobjective_hessian
ModelingToolkit.generate_control_jacobian
ModelingToolkit.generate_paramjac
ModelingToolkit.build_explicit_observed_function
ModelingToolkit.generate_control_function
ModelingToolkit.generate_update_A
ModelingToolkit.generate_update_b
ModelingToolkit.generate_semiquadratic_functions
ModelingToolkit.generate_semiquadratic_jacobian
ModelingToolkit.get_semiquadratic_W_sparsity
ModelingToolkitBase.CompilerOptions
ModelingToolkitBase.generate_custom_function
```

## Independent trailing arguments

[`ModelingToolkitBase.generate_custom_function`](@ref) accepts scalar or array symbolic
inputs in `extra_args`. These inputs are independent of the system's parameters and
follow the standard arguments, in tuple order. For a time-dependent system the resulting
functions take `(u, p, t, extra_args...)` and `(out, u, p, t, extra_args...)`; a
time-independent system omits `t`. The parameter argument retains the system's split
or unsplit representation.

```@example custom_function_extra_args
using ModelingToolkitBase

@independent_variables t
@variables x(t) dt a[1:2] b[1:2]
@parameters k = 2.0
sys = complete(System(Equation[], t, [x], [k]; name = :sys))
f, f! = generate_custom_function(
    sys, [x + dt * (k * a[1] + b[2]), x + a[2]];
    extra_args = (dt, a, b), expression = Val(false)
)
p = MTKParameters(sys, [x => 3.0])
u = [3.0]
a_values = [5.0, 7.0]
b_values = [11.0, 13.0]
result = f(u, p, 1.0, 0.5, a_values, b_values)
out = zeros(2)
f!(out, u, p, 1.0, 0.5, a_values, b_values)
@assert result == out == [14.5, 10.0]
result
```

With `expression = Val(true)`, the same argument order applies to the generated
expressions. The default `extra_args = ()` adds no arguments. Cache arguments supplied
through `cachesyms` retain their existing representation and position, including their
placement in the parameter argument for split systems. They are separate from `extra_args`.

For functions such as jacobian calculation which require symbolic computation, there
are `calculate_*` equivalents to obtain the symbolic result without building a function.

```@docs
ModelingToolkit.calculate_tgrad
ModelingToolkit.calculate_jacobian
ModelingToolkit.jacobian_dae_sparsity
ModelingToolkit.calculate_hessian
ModelingToolkit.calculate_massmatrix
ModelingToolkit.W_sparsity
ModelingToolkit.calculate_W_prototype
ModelingToolkit.calculate_cost_gradient
ModelingToolkit.calculate_cost_hessian
ModelingToolkit.cost_hessian_sparsity
ModelingToolkit.calculate_constraint_jacobian
ModelingToolkit.calculate_constraint_hessian
ModelingToolkit.calculate_multiobjective_jacobian
ModelingToolkit.calculate_multiobjective_hessian
ModelingToolkit.calculate_control_jacobian
ModelingToolkit.calculate_paramjac
ModelingToolkit.calculate_A_b
```

The constraint function and its derivatives are built from the residuals of the
constraints, with one row per element of each array-valued constraint.

```@docs
ModelingToolkitBase.canonical_constraints
ModelingToolkitBase.constraint_residual
```

A system can be marked as unsupported by symbolic automatic differentiation, in which case
the `calculate_*` functions above throw instead of producing a wrong derivative.

```@docs
ModelingToolkitBase.SymbolicADDisallowed
ModelingToolkitBase.check_symbolic_ad_allowed
```

All code generation eventually calls `build_function_wrapper`.

```@docs
ModelingToolkit.build_function_wrapper
```
