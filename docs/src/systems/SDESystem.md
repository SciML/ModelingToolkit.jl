# SDESystem

## System Constructors

```@docs
SDESystem
```

## Composition and Accessor Functions

- `sys.eqs` or `equations(sys)`: The equations that define the SDE.
- `sys.states` or `states(sys)`: The set of states in the SDE.
- `sys.parameters` or `parameters(sys)`: The parameters of the SDE.
- `sys.iv` or `independent_variable(sys)`: The independent variable of the SDE.

## Transformations

## Applicable Calculation and Generation Functions

```julia
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
SDEFunction
SDEProblem
```
