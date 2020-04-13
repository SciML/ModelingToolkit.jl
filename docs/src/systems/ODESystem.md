# ODESystem

## System Constructors

```@docs
ODESystem
```

## Composition and Accessor Functions

- `sys.eqs` or `equations(sys)`: The equations that define the ODE.
- `sys.states` or `states(sys)`: The set of states in the ODE.
- `sys.parameters` or `parameters(sys)`: The parameters of the ODE.
- `sys.iv` or `independent_variable(sys)`: The independent variable of the ODE.

## Transformations

```@docs
ode_order_lowering
```

## Function Calculation and Generation

```@docs
calculate_jacobian
calculate_tgrad
calculate_factorized_W
generate_jacobian
generate_tgrad
generate_factorized_W
ODEFunction
```

## Problem Constructors

```@docs
ODEProblem
```
