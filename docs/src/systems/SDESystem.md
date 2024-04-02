# SDESystem

## System Constructors

```@docs
SDESystem
```

To convert an `ODESystem` to an `SDESystem` directly:

```
ode = ODESystem(eqs,t,[x,y,z],[σ,ρ,β])
sde = SDESystem(ode, noiseeqs)
```

## Composition and Accessor Functions

  - `get_eqs(sys)` or `equations(sys)`: The equations that define the SDE.
  - `get_unknowns(sys)` or `unknowns(sys)`: The set of unknowns in the SDE.
  - `get_ps(sys)` or `parameters(sys)`: The parameters of the SDE.
  - `get_iv(sys)`: The independent variable of the SDE.
  - `continuous_events(sys)`: The set of continuous events in the SDE.
  - `discrete_events(sys)`: The set of discrete events in the SDE.
  - `alg_equations(sys)`: The algebraic equations (i.e. that does not contain a differential) that defines the ODE.
  - `get_alg_eqs(sys)`: The algebraic equations (i.e. that does not contain a differential) that defines the ODE. Only returns equations of the current-level system.
  - `diff_equations(sys)`: The differential equations (i.e. that contain a differential) that defines the ODE.
  - `get_diff_eqs(sys)`: The differential equations (i.e. that contain a differential) that defines the ODE. Only returns equations of the current-level system.
  - `has_alg_equations(sys)`: Returns `true` if the ODE contains any algebraic equations (i.e. that does not contain a differential).
  - `has_alg_eqs(sys)`: Returns `true` if the ODE contains any algebraic equations (i.e. that does not contain a differential). Only considers the current-level system.
  - `has_diff_equations(sys)`: Returns `true` if the ODE contains any differential equations (i.e. that does contain a differential).
  - `has_diff_eqs(sys)`: Returns `true` if the ODE contains any differential equations (i.e. that does contain a differential). Only considers the current-level system.

## Transformations

```@docs; canonical=false
structural_simplify
alias_elimination
```

```@docs
ModelingToolkit.Girsanov_transform
```

## Analyses

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

## Problem Constructors

```@docs
SDEFunction(sys::ModelingToolkit.SDESystem, args...)
SDEProblem(sys::ModelingToolkit.SDESystem, args...)
```

## Expression Constructors

```@docs
SDEFunctionExpr
SDEProblemExpr
```
