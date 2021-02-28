# ControlSystem

## System Constructors

```@docs
ControlSystem
```

## Composition and Accessor Functions

- `get_eqs(sys` or `equations(sys)`: The equations that define the system.
- `get_states(sys)` or `states(sys)`: The set of states in the system.
- `get_ps(sys)` or `parameters(sys)`: The parameters of the system.
- `get_controls(sys)` or `controls(sys)`: The control variables of the system

## Transformations

```@docs
ModelingToolkit.runge_kutta_discretize
structural_simplify
```
