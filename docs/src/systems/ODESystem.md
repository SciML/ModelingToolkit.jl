# ODESystem

## System Constructors

```@docs
ODESystem
```

## Composition and Accessor Functions

```@docs
canonical = false here
get_eqs
equations
get_unknowns
unknowns
get_ps
parameters
get_iv
get_u0_p
continuous_events
discrete_events
alg_equations
get_alg_eqs
diff_equations
get_diff_eqs
has_alg_equations
has_alg_eqs
has_diff_equations
has_diff_eqs


## Transformations

```@docs
structural_simplify
ode_order_lowering
dae_index_lowering
change_independent_variable
liouville_transform
alias_elimination
tearing
```

## Analyses

```@docs
ModelingToolkit.islinear
ModelingToolkit.isautonomous
ModelingToolkit.isaffine
```

## Applicable Calculation and Generation Functions

```@docs; canonical=false
calculate_jacobian
calculate_tgrad
calculate_factorized_W
generate_jacobian
generate_tgrad
generate_factorized_W
jacobian_sparsity
```

## Standard Problem Constructors

```@docs
ODEFunction(sys::ModelingToolkit.AbstractODESystem, args...)
ODEProblem(sys::ModelingToolkit.AbstractODESystem, args...)
SteadyStateProblem(sys::ModelingToolkit.AbstractODESystem, args...)
DAEProblem(sys::ModelingToolkit.AbstractODESystem, args...)
```

## Expression Constructors

```@docs
ODEFunctionExpr
DAEFunctionExpr
SteadyStateProblemExpr
```
