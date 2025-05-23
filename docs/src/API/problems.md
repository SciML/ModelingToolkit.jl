# Building and solving numerical problems

Systems are numerically solved by building and solving the appropriate problem type.
Numerical solvers expect to receive functions taking a predefeined set of arguments
and returning specific values. This format of argument and return value depends on
the function and the problem. ModelingToolkit is capable of compiling and generating
code for a variety of such numerical problems.

## Dynamical systems

```@docs
ODEFunction(::System, args...)
ODEProblem(::System, args...)
DAEFunction(::System, args...)
DAEProblem(::System, args...)
SDEFunction(::System, args...)
SDEProblem(::System, args...)
DDEFunction(::System, args...)
DDEProblem(::System, args...)
SDDEFunction(::System, args...)
SDDEProblem(::System, args...)
JumpProblem(::System, args...)
BVProblem(::System, args...)
DiscreteProblem(::System, args...)
ImplicitDiscreteProblem(::System, args...)
```

## Nonlinear systems

```@docs
NonlinearFunction(::System, args...)
NonlinearProblem(::System, args...)
SCCNonlinearProblem(::System, args...)
NonlinearLeastSquaresProblem(::System, args...)
SteadyStateProblem(::System, args...)
IntervalNonlinearFunction(::System, args...)
IntervalNonlinearProblem(::System, args...)
```

## Optimization and optimal control

```@docs
OptimizationFunction(::System, args...)
OptimizationProblem(::System, args...)
ODEInputFunction(::System, args...)
JuMPDynamicOptProblem(::System, args...)
InfiniteOptDynamicOptProblem,(::System, args...)
CasADiDynamicOptProblem(::System, args...)
DynamicOptSolution
```

## Linear analysis

```@docs
linearization_function
LinearizationProblem
linearize
CommonSolve.solve(::LinearizationProblem)
linearize_symbolic
```

There are also utilities for manipulating the results of these analyses in a symbolic context.

```@docs
ModelingToolkit.similarity_transform
ModelingToolkit.reorder_unknnowns
```

### Analysis point transformations

Linear analysis can also be done using analysis points to perform several common
workflows.

```@docs
get_sensitivity_function
get_sensitivity
get_comp_sensitivity_function
get_comp_sensitivity
get_looptransfer_function
get_looptransfer
open_loop
```
