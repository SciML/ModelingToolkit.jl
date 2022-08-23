# SDESystem

## System Constructors

```@docs
SDESystem
```

To convert an `ODESystem` to an `SDESystem` directly:
```
ode = ODESystem(eqs,t,[x,y,z],[σ,ρ,β])
sde = SDESystem(ode, noiseeqs)
```

## Composition and Accessor Functions

- `get_eqs(sys)` or `equations(sys)`: The equations that define the SDE.
- `get_states(sys)` or `states(sys)`: The set of states in the SDE.
- `get_ps(sys)` or `parameters(sys)`: The parameters of the SDE.
- `get_iv(sys)`: The independent variable of the SDE.

## Transformations

```@docs
structural_simplify
alias_elimination
Girsanov_transform
```

## Analyses

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
SDEFunction(sys::ModelingToolkit.SDESystem, args...)
SDEProblem(sys::ModelingToolkit.SDESystem, args...)
```

## Expression Constructors

```@docs
SDEFunctionExpr
SDEProblemExpr
```
