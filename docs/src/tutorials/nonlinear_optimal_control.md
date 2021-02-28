# Nonlinear Optimal Control

#### Note: this is still a work in progress!

The `ControlSystem` type is an interesting system because, unlike other
system types, it cannot be numerically solved on its own. Instead, it must be
transformed into another system before solving. Standard methods such as the
"direct method", "multiple shooting", or "discretize-then-optimize" can all be
phrased as symbolic transformations to a `ControlSystem`: this is the strategy
of this methodology.

## Defining a Nonlinear Optimal Control Problem

Here we will start by defining a classic optimal control problem. Let:

```math
x^{′′} = u^3(t)
```

where we want to optimize our controller `u(t)` such that the following is
minimized:

```math
L(\theta) = \sum_i \Vert 4 - x(t_i) \Vert + 2 \Vert x^\prime(t_i) \Vert + \Vert u(t_i) \Vert
```

where ``i`` is measured on (0,8) at 0.01 intervals. To do this, we rewrite the
ODE in first order form:

```math
\begin{aligned}
x^\prime &= v \\
v^′ &= u^3(t) \\
\end{aligned}
```

and thus

```math
L(\theta) = \sum_i \Vert 4 - x(t_i) \Vert + 2 \Vert v(t_i) \Vert + \Vert u(t_i) \Vert
```

is our loss function on the first order system.

Defining such a control system is similar to an `ODESystem`, except we must also
specify a control variable `u(t)` and a loss function. Together, this problem
looks as follows:

```julia
using ModelingToolkit

@variables t x(t) v(t) u(t)
@parameters p[1:2]
D = Differential(t)

loss = (4-x)^2 + 2v^2 + u^2
eqs = [
    D(x) ~ v - p[2]*x
    D(v) ~ p[1]*u^3 + v
]

sys = ControlSystem(loss,eqs,t,[x,v],[u],p)
```

## Solving a Control Problem via Discretize-Then-Optimize

One common way to solve nonlinear optimal control problems is by transforming
them into an optimization problem by performing a Runge-Kutta discretization
of the differential equation system and imposing equalities between variables
in the same steps. This can be done via the `runge_kutta_discretize` transformation
on the `ControlSystem`. While a tableau `tab` can be specified, it defaults to
a 5th order RadauIIA collocation, which is a common method in the field. To
perform this discretization, we simply need to give a `dt` and a timespan on which
to discretize:

```julia
dt = 0.1
tspan = (0.0,1.0)
sys = runge_kutta_discretize(sys,dt,tspan)
```

Now `sys` is an `OptimizationSystem` which, when solved, gives the values of
`x(t)`, `v(t)`, and `u(t)`. Thus we solve the `OptimizationSystem` using
GalacticOptim.jl:

```julia
u0 = rand(length(states(sys))) # guess for the state values
prob = OptimizationProblem(sys,u0,[0.1,0.1],grad=true)
sol = solve(prob,BFGS())
```

And this is missing some nice interfaces and ignores the equality constraints
right now so the tutorial is not complete.
