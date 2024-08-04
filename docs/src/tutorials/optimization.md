# Modeling Optimization Problems

ModelingToolkit.jl is not only useful for generating initial value problems (`ODEProblem`).
The package can also build optimization systems.

!!! note
    
    The high level `@mtkmodel` macro used in the
    [getting started tutorial](@ref getting_started)
    is not yet compatible with `OptimizationSystem`.
    We thus have to use a lower level interface to define optimization systems.
    For an introduction to this interface, read the
    [programmatically generating ODESystems tutorial](@ref programmatically).

## Unconstrained Rosenbrock Function

Let's optimize the classical _Rosenbrock function_ in two dimensions.

```@example optimization
using ModelingToolkit, Optimization, OptimizationOptimJL
@variables begin
    x = 1.0, [bounds = (-2.0, 2.0)]
    y = 3.0, [bounds = (-1.0, 3.0)]
end
@parameters a=1.0 b=1.0
rosenbrock = (a - x)^2 + b * (y - x^2)^2
@mtkbuild sys = OptimizationSystem(rosenbrock, [x, y], [a, b])
```

Every optimization problem consists of a set of optimization variables.
In this case, we create two variables: `x` and `y`,
with initial guesses `1` and `3` for their optimal values.
Additionally, we assign box constraints for each of them, using `bounds`,
Bounds is an example of symbolic metadata.
Fore more information, take a look at the symbolic metadata
[documentation page](@ref symbolic_metadata).

We also create two parameters with `@parameters`.
Parameters are useful if you want to solve the same optimization problem multiple times,
with different values for these parameters.
Default values for these parameters can also be assigned, here `1` is used for both `a` and `b`.
These optimization values and parameters are used in an objective function, here the Rosenbrock function.

Next, the actual `OptimizationProblem` can be created.
The initial guesses for the optimization variables can be overwritten, via an array of `Pairs`,
in the second argument of `OptimizationProblem`.
Values for the parameters of the system can also be overwritten from their default values,
in the third argument of `OptimizationProblem`.
ModelingToolkit is also capable of constructing analytical gradients and Hessians of the objective function.

```@example optimization
u0 = [y => 2.0]
p = [b => 100.0]

prob = OptimizationProblem(sys, u0, p, grad = true, hess = true)
u_opt = solve(prob, GradientDescent())
```

A visualization of the Rosenbrock function is depicted below.

```@example optimization
using Plots
x_plot = -2:0.01:2
y_plot = -1:0.01:3
contour(
    x_plot, y_plot, (x, y) -> (1 - x)^2 + 100 * (y - x^2)^2, fill = true, color = :viridis,
    ratio = :equal, xlims = (-2, 2))
scatter!([u_opt[1]], [u_opt[2]], ms = 10, label = "minimum")
```

## Rosenbrock Function with Constraints

ModelingToolkit is also capable of handing more complicated constraints than box constraints.
Non-linear equality and inequality constraints can be added to the `OptimizationSystem`.
Let's add an inequality constraint to the previous example:

```@example optimization_constrained
using ModelingToolkit, Optimization, OptimizationOptimJL

@variables begin
    x = 0.14, [bounds = (-2.0, 2.0)]
    y = 0.14, [bounds = (-1.0, 3.0)]
end
@parameters a=1.0 b=100.0
rosenbrock = (a - x)^2 + b * (y - x^2)^2
cons = [
    x^2 + y^2 ≲ 1
]
@mtkbuild sys = OptimizationSystem(rosenbrock, [x, y], [a, b], constraints = cons)
prob = OptimizationProblem(sys, [], grad = true, hess = true, cons_j = true, cons_h = true)
u_opt = solve(prob, IPNewton())
```

Inequality constraints are constructed via a `≲` (or `≳`).
[(To write these symbols in your own code write `\lesssim` or `\gtrsim` and then press tab.)]
(https://docs.julialang.org/en/v1/manual/unicode-input/)
An equality constraint can be specified via a `~`, e.g., `x^2 + y^2 ~ 1`.

A visualization of the Rosenbrock function and the inequality constraint is depicted below.

```@example optimization_constrained
using Plots
x_plot = -2:0.01:2
y_plot = -1:0.01:3
contour(
    x_plot, y_plot, (x, y) -> (1 - x)^2 + 100 * (y - x^2)^2, fill = true, color = :viridis,
    ratio = :equal, xlims = (-2, 2))
contour!(x_plot, y_plot, (x, y) -> x^2 + y^2, levels = [1], color = :lightblue, line = 4)
scatter!([u_opt[1]], [u_opt[2]], ms = 10, label = "minimum")
```

## Nested Systems

Needs more text, but it's super cool and auto-parallelizes and sparsifies too.
Plus, you can hierarchically nest systems to have it generate huge
optimization problems.
