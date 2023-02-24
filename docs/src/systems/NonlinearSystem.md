# NonlinearSystem

## System Constructors

```@docs
NonlinearSystem
```

## Composition and Accessor Functions

  - `get_eqs(sys)` or `equations(sys)`: The equations that define the nonlinear system.
  - `get_states(sys)` or `states(sys)`: The set of states in the nonlinear system.
  - `get_ps(sys)` or `parameters(sys)`: The parameters of the nonlinear system.
  - `get_u0_p(sys, u0map, parammap)` Numeric arrays for the initial condition and parameters given `var => value` maps.

## Transformations

```@docs
structural_simplify
alias_elimination
tearing
```

## Analyses

```@docs
ModelingToolkit.isaffine
ModelingToolkit.islinear
```

## Applicable Calculation and Generation Functions

```julia
calculate_jacobian
generate_jacobian
jacobian_sparsity
```

## Problem Constructors

```@docs
NonlinearFunction(sys::ModelingToolkit.NonlinearSystem, args...)
NonlinearProblem(sys::ModelingToolkit.NonlinearSystem, args...)
```

## Torn Problem Constructors

```@docs
BlockNonlinearProblem
```

## Expression Constructors

```@docs
NonlinearFunctionExpr
NonlinearProblemExpr
```
