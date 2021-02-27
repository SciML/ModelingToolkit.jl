# JumpSystem

## System Constructors

```@docs
JumpSystem
```

## Composition and Accessor Functions

- `get_eqs(sys)` or `equations(sys)`: The equations that define the jump system.
- `get_states(sys)` or `states(sys)`: The set of states in the jump system.
- `get_ps(sys)` or `parameters(sys)`: The parameters of the jump system.
- `independent_variable(sys)`: The independent variable of the jump system.

## Problem Constructors

```@docs
DiscreteProblem
JumpProblem
```
