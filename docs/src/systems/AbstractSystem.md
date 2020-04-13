# The AbstractSystem Interface

## Overview

The `AbstractSystem` interface is the core of the system level of ModelingToolkit.jl.
It establishes a common set of functionality that is used between systems
from ODEs and chemical reactions, allowing users to have a common framework for
model manipulation and compilation.

## Composition and Accessor Functions

Each `AbstractSystem` has lists of variables in context, such as distinguishing
parameters vs states. In addition, an `AbstractSystem` also can hold other
`AbstractSystem` types. Direct accessing of the values, such as `sys.states`,
gives the immediate list, while the accessor functions `states(sys)` gives the
total set which includes that of all systems held inside.

The values which are common to all `AbstractSystem`s are:

- `sys.eqs` or `equations(sys)`: The equations that define the system.
- `sys.states` or `states(sys)`: The set of states in the system.
- `sys.parameters` or `parameters(sys)`: The parameters of the system.
- `sys.systems`: The subsystems of the system.

## Transformations

Transformations are functions which send a valid `AbstractSystem` definition to
another `AbstractSystem`. These are passes like optimizations (ex: Block-Lower
Triangle transformations) or changes to the representation which allow for
alternative numerical methods to be utilized on the model (ex: DAE index reduction).

## Function Calculation and Generation

The calculation and generation functions allow for calculating additional
quantities to enhance the numerical methods applied to the resulting system.
The calculations, like `calculate_jacobian`, generate ModelingToolkit IR for
the Jacobian of the system, while the generations, like `generate_jacobian`,
generate compiled output for the numerical solvers by applying `build_function`
to the generated code. Additionally, many systems have function type outputs
which cobble together the generation functionality for a system, for example
`ODEFunction` can be used to generate a DifferentialEquations-based `ODEFunction`
with compiled version of the ODE itself, the Jacobian, the mass matrix, etc.

Below are the possible calculation and generation functions:

```@docs
calculate_tgrad
calculate_grad
calculate_jacobian
calculate_factorized_W
calculate_hessian
generate_jacobian
generate_grad
generate_tgrad
generate_factorized_W
generate_hessian
```

## Problem Constructors

At the end, the system types have `DEProblem` constructors, like `ODEProblem`,
which allow for directly generating the problem types required for numerical
methods. The first argument is always the `AbstractSystem`, and the proceeding
arguments match the argument order of their original constructors. Whenever an
array would normally be provided, such as `u0` the initial condition of an
`ODEProblem`, it is instead replaced with a variable map, i.e. an array of
pairs `var=>value` which allows the user to designate the values without having
to know the order that ModelingToolkit is internally using.
