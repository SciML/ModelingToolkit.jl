# Modeling Nonlinear Systems

ModelingToolkit.jl is not only useful for generating initial value problems (`ODEProblem`).
The package can also build nonlinear systems.
This is, for example, useful for finding the steady state of an ODE.
This steady state is reached when the nonlinear system of differential equations equals zero.

!!! note
    
    The high level `@mtkmodel` macro used in the
    [getting started tutorial](@ref getting_started)
    is not yet compatible with `NonlinearSystem`.
    We thus have to use a lower level interface to define nonlinear systems.
    For an introduction to this interface, read the
    [programmatically generating ODESystems tutorial](@ref programmatically).

```@example nonlinear
using ModelingToolkit, NonlinearSolve

# Define a nonlinear system
@variables x y z
@parameters σ ρ β
eqs = [0 ~ σ * (y - x)
       0 ~ x * (ρ - z) - y
       0 ~ x * y - β * z]
@mtkbuild ns = NonlinearSystem(eqs)

guesses = [x => 1.0, y => 0.0, z => 0.0]
ps = [σ => 10.0, ρ => 26.0, β => 8 / 3]

prob = NonlinearProblem(ns, guesses, ps)
sol = solve(prob, NewtonRaphson())
```

We found the `x`, `y` and `z` for which the right hand sides of `eqs` are all equal to zero.

Just like with `ODEProblem`s we can generate the `NonlinearProblem` with its analytical
Jacobian function:

```@example nonlinear
prob = NonlinearProblem(ns, guesses, ps, jac = true)
sol = solve(prob, NewtonRaphson())
```
