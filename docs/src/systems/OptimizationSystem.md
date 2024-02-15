# OptimizationSystem

## System Constructors

```@docs
OptimizationSystem
```

## Composition and Accessor Functions

  - `get_op(sys)`: The objective to be minimized.
  - `get_unknowns(sys)` or `unknowns(sys)`: The set of unknowns for the optimization.
  - `get_ps(sys)` or `parameters(sys)`: The parameters for the optimization.
  - `get_constraints(sys)` or `constraints(sys)`: The constraints for the optimization.

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
OptimizationProblem(sys::ModelingToolkit.OptimizationSystem, args...)
```

## Expression Constructors

```@docs
OptimizationProblemExpr
```
