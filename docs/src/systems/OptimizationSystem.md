# OptimizationSystem

## System Constructors

```@docs
OptimizationSystem
```

## Composition and Accessor Functions

- `sys.eqs` or `equations(sys)`: The equation to be minimized.
- `sys.states` or `states(sys)`: The set of states for the optimization.
- `sys.parameters` or `parameters(sys)`: The parameters for the optimization.

## Transformations

## Function Calculation and Generation

```@docs
calculate_grad
calculate_hes
generate_grad
generate_hes
```

## Problem Constructors

```@docs
OptimizationProblem
```
