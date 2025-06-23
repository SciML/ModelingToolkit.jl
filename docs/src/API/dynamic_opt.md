# [Dynamic Optimization Solvers](@id dynamic_opt_api)

Currently 4 backends are exposed for solving dynamic optimization problems using collocation: JuMP, InfiniteOpt, CasADi, and Pyomo.

Please note that there are differences in how to construct the collocation solver for the different cases. For example, the Python based ones, CasADi and Pyomo, expect the solver to be passed in as a string (CasADi and Pyomo come pre-loaded with Ipopt, but other solvers may need to be manually installed using `pip` or `conda`), while JuMP/InfiniteOpt expect the optimizer object to be passed in directly:

```
JuMPCollocation(Ipopt.Optimizer, constructRK4())
CasADiCollocation("ipopt", constructRK4())
```

**JuMP** and **CasADi** collocation require an ODE tableau to be passed in. These can be constructed by calling the `constructX()` functions from DiffEqDevTools. The list of tableaus can be found [here](https://docs.sciml.ai/DiffEqDevDocs/dev/internals/tableaus/). If none is passed in, both solvers will default to using Radau second-order with five collocation points.

**Pyomo** and **InfiniteOpt** each have their own built-in collocation methods.

 1. **InfiniteOpt**: The list of InfiniteOpt collocation methods can be found [in the table on this page](https://infiniteopt.github.io/InfiniteOpt.jl/stable/guide/derivative/). If none is passed in, the solver defaults to `FiniteDifference(Backward())`, which is effectively implicit Euler.
 2. **Pyomo**: The list of Pyomo collocation methods can be found [at the bottom of this page](https://github.com/SciML/Pyomo.jl). If none is passed in, the solver defaults to a `LagrangeRadau(3)`.

Some examples of the latter two collocations:

```julia
PyomoCollocation("ipopt", LagrangeRadau(2))
InfiniteOptCollocation(Ipopt.Optimizer, OrthogonalCollocation(3))
```

```@docs; canonical = false
JuMPCollocation
InfiniteOptCollocation
CasADiCollocation
PyomoCollocation
CommonSolve.solve(::AbstractDynamicOptProblem)
```

### Problem constructors

```@docs; canonical = false
JuMPDynamicOptProblem
InfiniteOptDynamicOptProblem
CasADiDynamicOptProblem
PyomoDynamicOptProblem
```
