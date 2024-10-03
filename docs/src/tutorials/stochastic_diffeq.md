# [Modeling with Stochasticity](@id SDE)

All models with `ODESystem` are deterministic. `SDESystem` adds another element
to the model: randomness. This is a
[stochastic differential equation](https://en.wikipedia.org/wiki/Stochastic_differential_equation)
which has a deterministic (drift) component and a stochastic (diffusion)
component. Let's take the Lorenz equation from the first tutorial and extend
it to have multiplicative noise by creating `@brownian` variables in the equations.

```@example SDE
using ModelingToolkit, StochasticDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D

# Define some variables
@parameters σ ρ β
@variables x(t) y(t) z(t)
@brownian a
eqs = [D(x) ~ σ * (y - x) + 0.1a * x,
    D(y) ~ x * (ρ - z) - y + 0.1a * y,
    D(z) ~ x * y - β * z + 0.1a * z]

@mtkbuild de = System(eqs, t)

u0map = [
    x => 1.0,
    y => 0.0,
    z => 0.0
]

parammap = [
    σ => 10.0,
    β => 26.0,
    ρ => 2.33
]

prob = SDEProblem(de, u0map, (0.0, 100.0), parammap)
sol = solve(prob, LambaEulerHeun())
```
