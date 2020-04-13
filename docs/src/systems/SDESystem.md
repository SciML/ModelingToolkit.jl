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

## Function Calculation and Generation

```@docs
calculate_jacobian(sys::ModelingToolkit.AbstractODESystem)
calculate_tgrad(sys::ModelingToolkit.AbstractODESystem)
calculate_factorized_W(sys::ModelingToolkit.AbstractODESystem, simplify)
generate_jacobian
generate_tgrad
generate_factorized_W
SDEFunction
```

## Problem Constructors

```@docs
SDEProblem
```
