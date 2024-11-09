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
where the magnitude of the noise scales with (0.3 times) the magnitude of each of the states:

```math
\begin{aligned}
\frac{dx}{dt} &= (\sigma (y-x))  &+ 0.3x\frac{dB}{dt} \\
\frac{dy}{dt} &= (x(\rho-z) - y) &+ 0.3y\frac{dB}{dt}  \\
\frac{dz}{dt} &= (xy - \beta z)  &+ 0.3z\frac{dB}{dt}  \\
\end{aligned}
```

Where $B$, is standard Brownian motion, also called the
[Wiener process](https://en.wikipedia.org/wiki/Wiener_process).
We use notation similar to the
[Langevin equation](https://en.wikipedia.org/wiki/Stochastic_differential_equation#Use_in_physics),
often used in physics.
By "multiplying" the equations by $dt$, the notation used in
[probability theory](https://en.wikipedia.org/wiki/Stochastic_differential_equation#Use_in_probability_and_mathematical_finance)
can be recovered.

We use this Langevin-like notation because it allows us to extend MTK modeling capacity from ODEs to SDEs,
using only a single new concept, `@brownian` variables, which represent $\frac{dB}{dt}$ in the above equation.

```@example SDE
using ModelingToolkit, StochasticDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
using Plots

@parameters σ=10.0 ρ=2.33 β=26.0
@variables x(t)=5.0 y(t)=5.0 z(t)=1.0
@brownian B
eqs = [D(x) ~ σ * (y - x) + 0.3x * B,
    D(y) ~ x * (ρ - z) - y + 0.3y * B,
    D(z) ~ x * y - β * z + 0.3z * B]

@mtkbuild de = System(eqs, t)
```

Even though we did not explicitly use `SDESystem`, ModelingToolkit can still infer this from the equations.

```@example SDE
typeof(de)
```

We continue by solving and plotting the SDE.

```@example SDE
prob = SDEProblem(de, [], (0.0, 100.0), [])
sol = solve(prob, SRIW1())
plot(sol, idxs = [(1, 2, 3)])
```

The noise present in all 3 equations is correlated, as can be seen on the below figure.
The figure also shows the multiplicative nature of the noise.
Because states `x` and `y` generally take on larger values,
the noise also takes on a more pronounced effect on these states compared to the state `z`.

```@example SDE
plot(sol)
```

If you want uncorrelated noise for each equation,
multiple `@brownian` variables have to be declared.

```@example SDE
@brownian Bx By Bz
eqs = [D(x) ~ σ * (y - x) + 0.3x * Bx,
    D(y) ~ x * (ρ - z) - y + 0.3y * By,
    D(z) ~ x * y - β * z + 0.3z * Bz]
@mtkbuild de = System(eqs, t)
prob = SDEProblem(de, [], (0.0, 100.0), [])
sol = solve(prob, SRIW1())
plot(sol)
```
