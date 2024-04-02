# ODESystem

## System Constructors

```@docs
ODESystem
```

## Composition and Accessor Functions

  - `get_eqs(sys)` or `equations(sys)`: The equations that define the ODE.
  - `get_unknowns(sys)` or `unknowns(sys)`: The set of unknowns in the ODE.
  - `get_ps(sys)` or `parameters(sys)`: The parameters of the ODE.
  - `get_iv(sys)`: The independent variable of the ODE.
  - `get_u0_p(sys, u0map, parammap)` Numeric arrays for the initial condition and parameters given `var => value` maps.
  - `continuous_events(sys)`: The set of continuous events in the ODE.
  - `discrete_events(sys)`: The set of discrete events in the ODE.
  - `alg_equations(sys)`: The algebraic equations (i.e. that does not contain a differential) that defines the ODE.
  - `get_alg_eqs(sys)`: The algebraic equations (i.e. that does not contain a differential) that defines the ODE. Only returns equations of the current-level system.
  - `diff_equations(sys)`: The differential equations (i.e. that contain a differential) that defines the ODE.
  - `get_diff_eqs(sys)`: The differential equations (i.e. that contain a differential) that defines the ODE. Only returns equations of the current-level system.
  - `has_alg_equations(sys)`: Returns `true` if the ODE contains any algebraic equations (i.e. that does not contain a differential).
  - `has_alg_eqs(sys)`: Returns `true` if the ODE contains any algebraic equations (i.e. that does not contain a differential). Only considers the current-level system.
  - `has_diff_equations(sys)`: Returns `true` if the ODE contains any differential equations (i.e. that does contain a differential).
  - `has_diff_eqs(sys)`: Returns `true` if the ODE contains any differential equations (i.e. that does contain a differential). Only considers the current-level system.

## Transformations

```@docs
structural_simplify
ode_order_lowering
dae_index_lowering
liouville_transform
alias_elimination
tearing
```

## Analyses

```@docs
ModelingToolkit.islinear
ModelingToolkit.isautonomous
ModelingToolkit.isaffine
```

## Applicable Calculation and Generation Functions

```@docs; canonical=false
calculate_jacobian
calculate_tgrad
calculate_factorized_W
generate_jacobian
generate_tgrad
generate_factorized_W
jacobian_sparsity
```

## Standard Problem Constructors

```@docs
ODEFunction(sys::ModelingToolkit.AbstractODESystem, args...)
ODEProblem(sys::ModelingToolkit.AbstractODESystem, args...)
SteadyStateProblem(sys::ModelingToolkit.AbstractODESystem, args...)
```

## Expression Constructors

```@docs
ODEFunctionExpr
DAEFunctionExpr
SteadyStateProblemExpr
```
