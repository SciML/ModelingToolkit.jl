# DiscreteSystem

## System Constructors

```@docs
DiscreteSystem
```

## Composition and Accessor Functions

  - `get_eqs(sys)` or `equations(sys)`: The equations that define the Discrete System.
  - `get_delay_val(sys)`: The delay of the Discrete System.
  - `get_iv(sys)`: The independent variable of the Discrete System.

## Transformations

## Analyses

## Applicable Calculation and Generation Functions

## Standard Problem Constructors

```@docs
DiscreteFunction(sys::ModelingToolkit.DiscreteSystem, args...)
DiscreteProblem(sys::ModelingToolkit.DiscreteSystem, args...)
```

## Expression Constructors

```@docs
DiscreteFunctionExpr
DiscreteProblemExpr
```
