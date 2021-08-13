# Modeling Nonlinear Systems

In this example we will go one step deeper and showcase the direct function
generation capabilities in ModelingToolkit.jl to build nonlinear systems.
Let's say we wanted to solve for the steady state of the previous ODE. This is
the nonlinear system defined by where the derivatives are zero. We use (unknown)
variables for our nonlinear system.

```julia
using ModelingToolkit, NonlinearSolve

@variables x y z
@parameters σ ρ β

# Define a nonlinear system
eqs = [0 ~ σ*(y-x),
       0 ~ x*(ρ-z)-y,
       0 ~ x*y - β*z]
@named ns = NonlinearSystem(eqs, [x,y,z], [σ,ρ,β])

guess = [x => 1.0,
         y => 0.0,
         z => 0.0]

ps = [
      σ => 10.0
      ρ => 26.0
      β => 8/3
      ]

prob = NonlinearProblem(ns,guess,ps)
sol = solve(prob,NewtonRaphson())
```

We can similarly ask to generate the `NonlinearProblem` with the analytical
Jacobian function:

```julia
prob = NonlinearProblem(ns,guess,ps,jac=true)
sol = solve(prob,NewtonRaphson())
```
