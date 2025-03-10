# ImplicitDiscreteSystem

## System Constructors

```@docs
ImplicitDiscreteSystem
```

## Composition and Accessor Functions

  - `get_eqs(sys)` or `equations(sys)`: The equations that define the implicit discrete system.
  - `get_unknowns(sys)` or `unknowns(sys)`: The set of unknowns in the implicit discrete system.
  - `get_ps(sys)` or `parameters(sys)`: The parameters of the implicit discrete system.
  - `get_iv(sys)`: The independent variable of the implicit discrete system
  - `discrete_events(sys)`: The set of discrete events in the implicit discrete system.

## Transformations

```@docs; canonical=false
structural_simplify
```

## Problem Constructors

```@docs; canonical=false
ImplicitDiscreteProblem(sys::ImplicitDiscreteSystem, u0map, tspan)
ImplicitDiscreteFunction(sys::ImplicitDiscreteSystem, args...)
```

## Discrete Domain

```@docs; canonical=false
Shift
```
