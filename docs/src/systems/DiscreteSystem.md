# DiscreteSystem

## System Constructors

```@docs
DiscreteSystem
```

## Composition and Accessor Functions

  - `get_eqs(sys)` or `equations(sys)`: The equations that define the discrete system.
  - `get_unknowns(sys)` or `unknowns(sys)`: The set of unknowns in the discrete system.
  - `get_ps(sys)` or `parameters(sys)`: The parameters of the discrete system.
  - `get_iv(sys)`: The independent variable of the discrete system
  - `discrete_events(sys)`: The set of discrete events in the discrete system.

## Transformations

```@docs; canonical=false
structural_simplify
```

## Problem Constructors

```@docs; canonical=false
DiscreteProblem(sys::DiscreteSystem, u0map, tspan)
DiscreteFunction(sys::DiscreteSystem, args...)
```
