# ControlSystem

## System Constructors

```@docs
ControlSystem
```

## Composition and Accessor Functions

- `sys.eqs` or `equations(sys)`: The equations that define the system.
- `sys.states` or `states(sys)`: The set of states in the system.
- `sys.parameters` or `parameters(sys)`: The parameters of the system.
- `sys.controls` or `controls(sys)`: The control variables of the system

## Transformations

```@docs
runge_kutta_discretize
```
