# Modeling with Stochasticity

All previous differential equations tutorials deal with deterministic `ODESystem`s.
In this tutorial, we add randomness.
In particular, we show how to represent a
[stochastic differential equation](https://en.wikipedia.org/wiki/Stochastic_differential_equation)
as a `SDESystem`.

!!! note
    
    The high level `@mtkmodel` macro used in the
    [getting started tutorial](@ref getting_started)
    is not yet compatible with `SDESystem`.
    We thus have to use a lower level interface to define stochastic differential equations.
    For an introduction to this interface, read the
    [programmatically generating ODESystems tutorial](@ref programmatically).

Let's take the Lorenz equation and add noise to each of the states.
To show the flexibility of ModelingToolkit,
we do not use homogeneous noise, with constant variance,
but instead use heterogeneous noise,
where the magnitude of the noise scales with (0.1 times) the magnitude of each of the states:

```math
\begin{aligned}
dx &= (\sigma (y-x))dt  &+ 0.1xdB \\
dy &= (x(\rho-z) - y)dt &+ 0.1ydB \\
dz &= (xy - \beta z)dt  &+ 0.1zdB \\
\end{aligned}
```

Where $B$, is standard Brownian motion, also called the
[Wiener process](https://en.wikipedia.org/wiki/Wiener_process).
In ModelingToolkit this can be done by `@brownian` variables.

```@example SDE
using ModelingToolkit, StochasticDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
using Plots

@parameters σ=10.0 ρ=2.33 β=26.0
@variables x(t)=5.0 y(t)=5.0 z(t)=1.0
@brownian B
eqs = [D(x) ~ σ * (y - x) + 0.1B * x,
    D(y) ~ x * (ρ - z) - y + 0.1B * y,
    D(z) ~ x * y - β * z + 0.1B * z]

@mtkbuild de = System(eqs, t)
```

Even though we did not explicitly use `SDESystem`, ModelingToolkit can still infer this from the equations.

```@example SDE
typeof(de)
```

We continue by solving and plotting the SDE.

```@example SDE
prob = SDEProblem(de, [], (0.0, 100.0), [])
sol = solve(prob, LambaEulerHeun())
plot(sol, idxs = [(1, 2, 3)])
```

The noise present in all 3 equations is correlated, as can be seen on the below figure.
If you want uncorrelated noise for each equation,
multiple `@brownian` variables have to be declared.

```@example SDE
@brownian Bx By Bz
```

The figure also shows the multiplicative nature of the noise.
Because states `x` and `y` generally take on larger values,
the noise also takes on a more pronounced effect on these states compared to the state `z`.

```@example SDE
plot(sol)
```
