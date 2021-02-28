# ODESystem

## System Constructors

```@docs
ODESystem
```

## Composition and Accessor Functions

- `get_eqs(sys)` or `equations(sys)`: The equations that define the ODE.
- `get_states(sys)` or `states(sys)`: The set of states in the ODE.
- `get_ps(sys)` or `parameters(sys)`: The parameters of the ODE.
- `independent_variable(sys)`: The independent variable of the ODE.

## Transformations

```@docs
structural_simplify
ode_order_lowering
liouville_transform
```

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
ODEFunction
ODEProblem
```
