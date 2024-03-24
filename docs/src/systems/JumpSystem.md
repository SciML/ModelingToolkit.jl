# JumpSystem

## System Constructors

```@docs
JumpSystem
```

## Composition and Accessor Functions

  - `get_eqs(sys)` or `equations(sys)`: The equations that define the jump system.
  - `get_unknowns(sys)` or `unknowns(sys)`: The set of unknowns in the jump system.
  - `get_ps(sys)` or `parameters(sys)`: The parameters of the jump system.
  - `get_iv(sys)`: The independent variable of the jump system.
  - `discrete_events(sys)`: The set of discrete events in the jump system.

## Transformations

```@docs; canonical=false
structural_simplify
```

## Analyses

## Problem Constructors

```@docs; canonical=false
DiscreteProblem(sys::JumpSystem, u0map, tspan)
```

```@docs
JumpProblem(sys::JumpSystem, prob, aggregator)
```
