# [Symbolic-Numeric Perturbation Theory for ODEs](@id perturb_diff)

In the [Mixed Symbolic-Numeric Perturbation Theory tutorial](https://symbolics.juliasymbolics.org/stable/tutorials/perturbation/), we discussed how to solve algebraic equations using **Symbolics.jl**. Here we extend the method to differential equations. The procedure is similar, but the Taylor series coefficients now become functions of an independent variable (usually time).

## Free fall in a varying gravitational field

Our first ODE example is a well-known physics problem: what is the altitude $x(t)$ of an object (say, a ball or a rocket) thrown vertically with initial velocity $ẋ(0)$ from the surface of a planet with mass $M$ and radius $R$? According to Newton's second law and law of gravity, it is the solution of the ODE

```math
ẍ = -\frac{GM}{(R+x)^2} = -\frac{GM}{R^2} \frac{1}{\left(1+ϵ\frac{x}{R}\right)^2}.
```

In the last equality, we introduced a perturbative expansion parameter $ϵ$. When $ϵ=1$, we recover the original problem. When $ϵ=0$, the problem reduces to the trivial problem $ẍ = -g$ with constant gravitational acceleration $g = GM/R^2$ and solution $x(t) = x(0) + ẋ(0) t - \frac{1}{2} g t^2$. This is a good setup for perturbation theory.

To make the problem dimensionless, we redefine $x \leftarrow x / R$ and $t \leftarrow t / \sqrt{R^3/GM}$. Then the ODE becomes

```@example perturbation
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
@variables ϵ x(t)
eq = D(D(x)) ~ -(1 + ϵ * x)^(-2)
```

Next, expand $x(t)$ in a series up to second order in $ϵ$:

```@example perturbation
using Symbolics
@variables y(t)[0:2] # coefficients
x_series = series(y, ϵ)
```

Insert this into the equation and collect perturbed equations to each order:

```@example perturbation
eq_pert = substitute(eq, x => x_series)
eqs_pert = taylor_coeff(eq_pert, ϵ, 0:2)
```

!!! note
    
    The 0-th order equation can be solved analytically, but ModelingToolkit does currently not feature automatic analytical solution of ODEs, so we proceed with solving it numerically.

These are the ODEs we want to solve. Now construct an `ODESystem`, which automatically inserts dummy derivatives for the velocities:

```@example perturbation
@mtkbuild sys = ODESystem(eqs_pert, t)
```

To solve the `ODESystem`, we generate an `ODEProblem` with initial conditions $x(0) = 0$, and $ẋ(0) = 1$, and solve it:

```@example perturbation
using OrdinaryDiffEq
u0 = Dict([unknowns(sys) .=> 0.0; D(y[0]) => 1.0]) # nonzero initial velocity
prob = ODEProblem(sys, u0, (0.0, 3.0))
sol = solve(prob)
```

This is the solution for the coefficients in the series for $x(t)$ and their derivatives. Finally, we calculate the solution to the original problem by summing the series for different $ϵ$:

```@example perturbation
using Plots
p = plot()
for ϵᵢ in 0.0:0.1:1.0
    plot!(p, sol, idxs = substitute(x_series, ϵ => ϵᵢ), label = "ϵ = $ϵᵢ")
end
p
```

This makes sense: for larger $ϵ$, gravity weakens with altitude, and the trajectory goes higher for a fixed initial velocity.

An advantage of the perturbative method is that we run the ODE solver only once and calculate trajectories for several $ϵ$ for free. Had we solved the full unperturbed ODE directly, we would need to do repeat it for every $ϵ$.

## Weakly nonlinear oscillator

Our second example applies perturbation theory to nonlinear oscillators -- a very important class of problems. As we will see, perturbation theory has difficulty providing a good solution to this problem, but the process is nevertheless instructive. This example closely follows chapter 7.6 of *Nonlinear Dynamics and Chaos* by Steven Strogatz.

The goal is to solve the ODE

```@example perturbation
eq = D(D(x)) + 2 * ϵ * D(x) + x ~ 0
```

with initial conditions $x(0) = 0$ and $ẋ(0) = 1$. With $ϵ = 0$, the problem reduces to the simple linear harmonic oscillator with the exact solution $x(t) = \sin(t)$.

We follow the same steps as in the previous example to construct the `ODESystem`:

```@example perturbation
eq_pert = substitute(eq, x => x_series)
eqs_pert = taylor_coeff(eq_pert, ϵ, 0:2)
@mtkbuild sys = ODESystem(eqs_pert, t)
```

We solve and plot it as in the previous example, and compare the solution with $ϵ=0.1$ to the exact solution $x(t, ϵ) = e^{-ϵ t} \sin(\sqrt{(1-ϵ^2)}\,t) / \sqrt{1-ϵ^2}$ of the unperturbed equation:

```@example perturbation
u0 = Dict([unknowns(sys) .=> 0.0; D(y[0]) => 1.0]) # nonzero initial velocity
prob = ODEProblem(sys, u0, (0.0, 50.0))
sol = solve(prob)
plot(sol, idxs = substitute(x_series, ϵ => 0.1); label = "Perturbative (ϵ=0.1)")

x_exact(t, ϵ) = exp(-ϵ * t) * sin(√(1 - ϵ^2) * t) / √(1 - ϵ^2)
plot!(sol.t, x_exact.(sol.t, 0.1); label = "Exact (ϵ=0.1)")
```

This is similar to Figure 7.6.2 in *Nonlinear Dynamics and Chaos*. The two curves fit well for the first couple of cycles, but then the perturbative solution diverges from the exact solution. The main reason is that the problem has two or more time-scales that introduce secular terms in the solution. One solution is to explicitly account for the two time scales and use an analytic method called *two-timing*, but this is outside the scope of this example.
