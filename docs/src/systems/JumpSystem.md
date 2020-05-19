# JumpSystem

## System Constructors

```@docs
JumpSystem
```

## Composition and Accessor Functions

- `sys.eqs` or `equations(sys)`: The equations that define the jump system.
- `sys.states` or `states(sys)`: The set of states in the jump system.
- `sys.parameters` or `parameters(sys)`: The parameters of the jump system.
- `sys.iv` or `independent_variable(sys)`: The independent variable of the jump system.

## Problem Constructors

```@docs
DiscreteProblem
JumpProblem
```
