# Solving Dynamic Optimization Problems
Systems in ModelingToolkit.jl can be directly converted to dynamic optimization or optimal control problems. In such systems, one has one or more input variables that are externally controlled to control the dynamics of the system. A dynamic optimization solves for the optimal time trajectory of the input variables in order to maximize or minimize a desired objective function. For example, a car driver might like to know how to step on the accelerator if the goal is to finish a race while using the least gas.

To begin, let us take a rocket launch example. The input variable here is the thrust exerted by the engine. The rocket state is described by its current height and velocity.
```julia
using ModelingToolkit
t = ModelingToolkit.t_nounits
D = ModelingToolkit.D_nounits

@parameters h_c m₀ h₀ g₀ D_c c Tₘ m_c
@variables begin
    h(..) 
    v(..) 
    m(..) [bounds = (m_c, 1)] 
    T(..) [input = true, bounds = (0, Tₘ)]
end

drag(h, v) = D_c * v^2 * exp(-h_c * (h - h₀) / h₀)
gravity(h) = g₀ * (h₀ / h)

eqs = [D(h(t)) ~ v(t),
       D(v(t)) ~ (T(t) - drag(h(t), v(t))) / m(t) - gravity(h(t)),
       D(m(t)) ~ -T(t) / c]

(ts, te) = (0.0, 0.2)
costs = [-h(te)]
cons = [T(te) ~ 0, m(te) ~ m_c]

@named rocket = ODESystem(eqs, t; costs, constraints = cons)
rocket, input_idxs = structural_simplify(rocket, ([T(t)], []))

u0map = [h(t) => h₀, m(t) => m₀, v(t) => 0]
pmap = [
    g₀ => 1, m₀ => 1.0, h_c => 500, c => 0.5 * √(g₀ * h₀), D_c => 0.5 * 620 * m₀ / g₀,
    Tₘ => 3.5 * g₀ * m₀, T(t) => 0.0, h₀ => 1, m_c => 0.6]
```
What we would like to optimize here is the final height of the rocket. We do this by providing a vector of expressions corresponding to the costs. By default, the sense of the optimization is to minimize the provided cost. So to maximize the rocket height at the final time, we write `-h(te)` as the cost.

Now we can construct a problem and solve it. Let us use JuMP as our backend here.
```julia
jprob = JuMPDynamicOptProblem(rocket, u0map, (ts, te), pmap; dt = 0.001, cse = false)
jsol = solve(jprob, JuMPCollocation(Ipopt.Optimizer, constructRadauIIA5()))
```

Let's plot our final solution and the controller here:
```julia
```

###### Free final time problems
There are additionally a class of dynamic optimization problems where we would like to know how to control our system to achieve something in the least time. Such problems are called free final time problems, since the final time is unknown. To model these problems in ModelingToolkit, we declare the final time as a parameter.

```julia
@variables x(..) v(..)
@variables u(..) [bounds = (-1.0, 1.0), input = true]
@parameters tf

constr = [v(tf) ~ 0, x(tf) ~ 0]
cost = [tf] # Minimize time

@named block = ODESystem(
    [D(x(t)) ~ v(t), D(v(t)) ~ u(t)], t; costs = cost, constraints = constr)

block, input_idxs = structural_simplify(block, ([u(t)], []))

u0map = [x(t) => 1.0, v(t) => 0.0]
tspan = (0.0, tf)
parammap = [u(t) => 0.0, tf => 1.0]
```

Please note that, at the moment, free final time problems cannot support constraints defined at definite time values, like `x(3) ~ 2`.

Let's plot our final solution and the controller for the block:
```julia
```

### Solvers
Currently 4 backends are exposed for solving dynamic optimization problems using collocation: JuMP, InfiniteOpt, CasADi, and Pyomo.

Please note that there are differences in how to construct the collocation solver for the different cases. For example, the Python based ones, CasADi and Pyomo, expect the solver to be passed in as a string (CasADi and Pyomo come pre-loaded with certain solvers, but other solvers may need to be manually installed using `pip` or `conda`), while JuMP/InfiniteOpt expect the optimizer object to be passed in directly:
```
JuMPCollocation(Ipopt.Optimizer, constructRK4())
CasADiCollocation("ipopt", constructRK4())
```

**JuMP** and **CasADi** collocation require an ODE tableau to be passed in. These can be constructed by calling the `constructX()` functions from DiffEqDevTools. If none is passed in, both solvers will default to using Radau second-order with five collocation points.

**Pyomo** and **InfiniteOpt** each have their own built-in collocation methods.
1. **InfiniteOpt**: The list of InfiniteOpt collocation methods can be found [in the table on this page](https://infiniteopt.github.io/InfiniteOpt.jl/stable/guide/derivative/). If none is passed in, the solver defaults to `FiniteDifference(Backward())`, which is effectively implicit Euler.
2. **Pyomo**: The list of Pyomo collocation methods can be found [here](). If none is passed in, the solver defaults to a `LagrangeRadau(3)`.

```@docs; canonical = false
JuMPCollocation
InfiniteOptCollocation
CasADiCollocation
PyomoCollocation
solve(::AbstractDynamicOptProblem)
```

### Problem constructors
```@docs; canonical = false
JuMPDynamicOptProblem
InfiniteOptDynamicOptProblem
CasADiDynamicOptProblem
PyomoDynamicOptProblem
```
