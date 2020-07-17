# NonlinearSystem

## System Constructors

```@docs
NonlinearSystem
```

## Composition and Accessor Functions

- `sys.eqs` or `equations(sys)`: The equations that define the nonlinear system.
- `sys.states` or `states(sys)`: The set of states in the nonlinear system.
- `sys.parameters` or `parameters(sys)`: The parameters of the nonlinear system.

## Transformations

## Applicable Calculation and Generation Functions

```julia
calculate_jacobian
generate_jacobian
jacobian_sparsity
```

## Problem Constructors

```@docs
NonlinearProblem
```
