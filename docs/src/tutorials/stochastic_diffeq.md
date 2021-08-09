# Modeling with Stochasticity

All models with `ODESystem` are deterministic. `SDESystem` adds another element
to the model: randomness. This is a
[stochastic differential equation](https://en.wikipedia.org/wiki/Stochastic_differential_equation)
which has a deterministic (drift) component and a stochastic (diffusion)
component. Let's take the Lorenz equation from the first tutorial and extend
it to have multiplicative noise.

```julia
using ModelingToolkit, StochasticDiffEq

# Define some variables
@parameters σ ρ β
@variables t x(t) y(t) z(t)
D = Differential(t)

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]

noiseeqs = [0.1*x,
            0.1*y,
            0.1*z]

de = SDESystem(eqs,noiseeqs,t,[x,y,z],[σ,ρ,β])

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

prob = SDEProblem(de,u0map,(0.0,100.0),parammap)
sol = solve(prob,SOSRI())
```
