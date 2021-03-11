# OptimizationSystem

## System Constructors

```@docs
OptimizationSystem
```

## Composition and Accessor Functions

- `get_eqs(sys)` or `equations(sys)`: The equation to be minimized.
- `get_states(sys)` or `states(sys)`: The set of states for the optimization.
- `get_ps(sys)` or `parameters(sys)`: The parameters for the optimization.

## Transformations

## Analyses

## Applicable Calculation and Generation Functions

```julia
calculate_gradient
calculate_hessian
generate_gradient
generate_hessian
hessian_sparsity
```

## Problem Constructors

```@docs
OptimizationProblem
```
