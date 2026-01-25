# [Symbolic-Numeric Perturbation Theory for ODEs](@id perturb_diff)

In the [Mixed Symbolic-Numeric Perturbation Theory tutorial](https://symbolics.juliasymbolics.org/stable/tutorials/perturbation/), we discussed how to solve algebraic equations using **Symbolics.jl**. Here we extend the method to differential equations. The procedure is similar, but the Taylor series coefficients now become functions of an independent variable (usually time).

## Free fall in a varying gravitational field

Our first ODE example is a well-known physics problem: what is the altitude $x(t)$ of an object (say, a ball or a rocket) thrown vertically with initial velocity $ẋ(0)$ from the surface of a planet with mass $M$ and radius $R$? According to Newton's second law and law of gravity, it is the solution of the ODE

```math
ẍ = -\frac{GM}{(R+x)^2} = -\frac{GM}{R^2} \frac{1}{\left(1+ϵ\frac{x}{R}\right)^2}.
```

In the last equality, we introduced a perturbative expansion parameter $ϵ$. When $ϵ=1$, we recover the original problem. When $ϵ=0$, the problem reduces to the trivial problem $ẍ = -g$ with constant gravitational acceleration $g = GM/R^2$ and solution $x(t) = x(0) + ẋ(0) t - \frac{1}{2} g t^2$. This is a good setup for perturbation theory.

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
@variables y₀(t) y₁(t) y₂(t) # coefficients (y₀=0th order, y₁=1st, y₂=2nd)
x_series = series([y₀, y₁, y₂], ϵ)
```

Insert this into the equation and collect perturbed equations to each order:

```@example perturbation
eq_pert = substitute(eq, x => x_series)
eqs_pert = taylor_coeff(eq_pert, ϵ, 0:2)
```

!!! note

    The 0-th order equation can be solved analytically, but ModelingToolkit does currently not feature automatic analytical solution of ODEs, so we proceed with solving it numerically.

These are the ODEs we want to solve. Now construct an `System`, which automatically inserts dummy derivatives for the velocities:

```@example perturbation
@mtkcompile sys = System(eqs_pert, t)
```

To solve the `System`, we generate an `ODEProblem` with initial conditions $x(0) = 0$, and $ẋ(0) = 1$, and solve it:

```@example perturbation
using OrdinaryDiffEq
# Get the system's unknowns for proper symbolic lookup
unks = unknowns(sys)
y₀_sys = first(filter(u -> contains(string(u), "y₀(t)"), unks))
y₁_sys = first(filter(u -> contains(string(u), "y₁(t)"), unks))
y₂_sys = first(filter(u -> contains(string(u), "y₂(t)"), unks))
u0 = Dict([unknowns(sys) .=> 0.0; D(y₀_sys) => 1.0]) # nonzero initial velocity
prob = ODEProblem(sys, u0, (0.0, 3.0))
sol = solve(prob)
```

This is the solution for the coefficients in the series for $x(t)$ and their derivatives. Finally, we calculate the solution to the original problem by summing the series for different $ϵ$:

```@example perturbation
using Plots
p = plot()
for ϵᵢ in 0.0:0.1:1.0
    # Compute x(t) = y₀(t) + y₁(t)*ϵ + y₂(t)*ϵ² from solution values
    xs = sol[y₀_sys] .+ ϵᵢ .* sol[y₁_sys] .+ ϵᵢ^2 .* sol[y₂_sys]
    plot!(p, sol.t, xs, label = "ϵ = $ϵᵢ")
end
p
```

This makes sense: for larger $ϵ$, gravity weakens with altitude, and the trajectory goes higher for a fixed initial velocity.

An advantage of the perturbative method is that we run the ODE solver only once and calculate trajectories for several $ϵ$ for free. Had we solved the full unperturbed ODE directly, we would need to do repeat it for every $ϵ$.

## Weakly nonlinear oscillator

Our second example applies perturbation theory to nonlinear oscillators -- a very important class of problems. As we will see, perturbation theory has difficulty providing a good solution to this problem, but the process is nevertheless instructive. This example closely follows chapter 7.6 of *Nonlinear Dynamics and Chaos* by Steven Strogatz.

The goal is to solve the ODE

```@example perturbation
eq2 = D(D(x)) + 2 * ϵ * D(x) + x ~ 0
```

with initial conditions $x(0) = 0$ and $ẋ(0) = 1$. With $ϵ = 0$, the problem reduces to the simple linear harmonic oscillator with the exact solution $x(t) = \sin(t)$.

We follow the same steps as in the previous example to construct the `System`, using new coefficients `z₀`, `z₁`, `z₂`:

```@example perturbation
@variables z₀(t) z₁(t) z₂(t) # new coefficients for this example
x_series2 = series([z₀, z₁, z₂], ϵ)
eq2_pert = substitute(eq2, x => x_series2)
eqs2_pert = taylor_coeff(eq2_pert, ϵ, 0:2)
@mtkcompile sys2 = System(eqs2_pert, t)

# Get the system's unknowns for proper symbolic lookup
unks2 = unknowns(sys2)
z₀_sys = first(filter(u -> contains(string(u), "z₀(t)"), unks2))
z₁_sys = first(filter(u -> contains(string(u), "z₁(t)"), unks2))
z₂_sys = first(filter(u -> contains(string(u), "z₂(t)"), unks2))
```

We solve and plot it as in the previous example, and compare the solution with $ϵ=0.1$ to the exact solution $x(t, ϵ) = e^{-ϵ t} \sin(\sqrt{(1-ϵ^2)}\,t) / \sqrt{1-ϵ^2}$ of the unperturbed equation:

```@example perturbation
u0 = [z₀_sys => 0.0, z₁_sys => 0.0, z₂_sys => 0.0, D(z₀_sys) => 1.0, D(z₁_sys) => 0.0, D(z₂_sys) => 0.0] # nonzero initial velocity
prob = ODEProblem(sys2, u0, (0.0, 50.0))
sol = solve(prob)
# Compute x(t) = z₀(t) + z₁(t)*ϵ + z₂(t)*ϵ² at ϵ=0.1
x_pert = sol[z₀_sys] .+ 0.1 .* sol[z₁_sys] .+ 0.1^2 .* sol[z₂_sys]
plot(sol.t, x_pert; label = "Perturbative (ϵ=0.1)")

x_exact(t, ϵ) = exp(-ϵ * t) * sin(√(1 - ϵ^2) * t) / √(1 - ϵ^2)
x_pert_at_pi2 = sol(π/2)[z₀_sys] + 0.1 * sol(π/2)[z₁_sys] + 0.1^2 * sol(π/2)[z₂_sys]
@assert isapprox(x_pert_at_pi2, x_exact(π/2, 0.1); atol = 1e-2) # compare around 1st peak # hide
plot!(sol.t, x_exact.(sol.t, 0.1); label = "Exact (ϵ=0.1)")
```

This is similar to Figure 7.6.2 in *Nonlinear Dynamics and Chaos*. The two curves fit well for the first couple of cycles, but then the perturbative solution diverges from the exact solution. The main reason is that the problem has two or more time-scales that introduce secular terms in the solution. One solution is to explicitly account for the two time scales and use an analytic method called *two-timing*, but this is outside the scope of this example.
