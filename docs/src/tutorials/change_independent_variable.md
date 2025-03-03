# Changing the independent variable of ODEs

Ordinary differential equations describe the rate of change of some dependent variables with respect to one independent variable.
For the modeler it is often most natural to write down the equations with a particular independent variable, say time $t$.
In many cases there are good reasons for reparametrizing ODEs in terms of a different independent variable:

1. One may want $y(x)$ as a function of $x$ instead of $(x(t), y(t))$ as a function of $t$
2. Some differential equations vary more nicely (e.g. less stiff or better behaved) with respect to one independent variable than another.
3. It can reduce the number of equations that must be solved (e.g. $y(x)$ is one equation, while $(x(t), y(t))$ are two).

To manually change the independent variable of an ODE, one must rewrite all equations in terms of a new variable and transform differentials with the chain rule.
This is mechanical and error-prone.
ModelingToolkit provides the utility function [`change_independent_variable`](@ref) that automates this process.

## 1. Get one dependent variable as a function of another

Consider a projectile shot with some initial velocity in a gravitational field.
```@example changeivar
using ModelingToolkit
@independent_variables t
D = Differential(t)
@variables x(t) y(t)
@parameters g = 9.81 v # gravitational acceleration and constant horizontal velocity
M1 = ODESystem([
    D(D(y)) ~ -g, D(x) ~ v # constant horizontal velocity
], t; defaults = [
    y => 0.0
], initialization_eqs = [
    #x ~ 0.0, # TODO: handle?
    D(x) ~ D(y) # equal initial horizontal and vertical velocity (45 °)
], name = :M) |> complete
M1s = structural_simplify(M1)
```
This is the standard parametrization that arises naturally from kinematics and Newton's laws.
It expresses the position $(x(t), y(t))$ as a function of time $t$.
But suppose we want to determine whether the projectile hits a target 10 meters away.
There are at least three ways of answering this:
* Solve the ODE for $(x(t), y(t))$ and use a callback to terminate when $x$ reaches 10 meters, and evaluate $y$ at the final time.
* Solve the ODE for $(x(t), y(t))$ and use root finding to find the time when $x$ reaches 10 meters, and evaluate $y$ at that time.
* Solve the ODE for $y(x)$ and evaluate it at 10 meters.

We will demonstrate the last method by changing the independent variable from $t$ to $x$.
This transformation is well-defined for any non-zero horizontal velocity $v$.
```@example changeivar
M2 = change_independent_variable(M1, M1.x; dummies = true) |> complete
@assert M2.x isa Num # hide
M2s = structural_simplify(M2; allow_symbolic = true)
```
The derivatives are now with respect to the new independent variable $x$, which can be accessed with `M2.x`.
Notice how the number of equations has also decreased from three to two, as $\mathrm{d}x/\mathrm{d}t$ has been turned into an observed equation.
It is straightforward to evolve the ODE for 10 meters and plot the resulting trajectory $y(x)$:
```@example changeivar
using OrdinaryDiffEq, Plots
prob = ODEProblem(M2s, [], [0.0, 10.0], [v => 8.0]) # throw 10 meters with x-velocity 8 m/s
sol = solve(prob, Tsit5())
plot(sol; idxs = M2.y) # must index by M2.y = y(x); not M1.y = y(t)!
```

## 2. Reduce stiffness by changing to a logarithmic time axis

In cosmology, the [Friedmann equations](https://en.wikipedia.org/wiki/Friedmann_equations) describe the expansion of the universe.
In terms of conformal time $t$, they can be written
```@example changeivar
@variables a(t) Ω(t) Ωr(t) Ωm(t) ΩΛ(t)
M1 = ODESystem([
    D(a) ~ √(Ω) * a^2
    Ω ~ Ωr + Ωm + ΩΛ
    D(Ωm) ~ -3 * D(a)/a * Ωm
    D(Ωr) ~ -4 * D(a)/a * Ωr
    D(ΩΛ) ~ 0
], t; initialization_eqs = [
    ΩΛ + Ωr + Ωm ~ 1
], name = :M) |> complete
M1s = structural_simplify(M1)
```
Of course, we can solve this ODE as it is:
```@example changeivar
prob = ODEProblem(M1s, [M1.a => 1.0, M1.Ωr => 5e-5, M1.Ωm => 0.3], (0.0, -3.3), [])
sol = solve(prob, Tsit5(); reltol = 1e-5)
@assert Symbol(sol.retcode) == :Unstable # surrounding text assumes this was unstable # hide
plot(sol, idxs = [M1.a, M1.Ωr/M1.Ω, M1.Ωm/M1.Ω, M1.ΩΛ/M1.Ω])
```
The solver becomes unstable due to stiffness.
Also notice the interesting dynamics taking place towards the end of the integration (in the early universe), which gets compressed into a very small time interval.
These ODEs would benefit from being defined with respect to a logarithmic "time" that better captures the evolution of the universe through *orders of magnitude* of time.

It is therefore common to write these ODEs in terms of $b = \ln a$.
To do this, we will change the independent variable in two stages; from $t$ to $a$ to $b$.
Notice that $\mathrm{d}a/\mathrm{d}t > 0$ provided that $\Omega > 0$, and $\mathrm{d}b/\mathrm{d}a > 0$, so the transformation is well-defined.
First, we transform from $t$ to $a$:
```@example changeivar
M2 = change_independent_variable(M1, M1.a; dummies = true) |> complete
```
Unlike the original, notice that this system is *non-autonomous* because the independent variable $a$ appears explicitly in the equations!
This means that to change the independent variable from $a$ to $b$, we must provide not only the rate of change relation $db(a)/da = \exp(-b)$, but *also* the equation $a(b) = \exp(b)$ so $a$ can be eliminated in favor of $b$:
```@example changeivar
a = M2.a
Da = Differential(a)
@variables b(a)
M3 = change_independent_variable(M2, b, [Da(b) ~ exp(-b), a ~ exp(b)]) |> complete
```
We can now solve and plot the ODE in terms of $b$:
```@example changeivar
M3s = structural_simplify(M3; allow_symbolic = true)
prob = ODEProblem(M3s, [M3.Ωr => 5e-5, M3.Ωm => 0.3], (0, -15), [])
sol = solve(prob, Tsit5())
@assert Symbol(sol.retcode) == :Success # surrounding text assumes the solution was successful # hide
plot(sol, idxs = [M3.Ωr/M3.Ω, M3.Ωm/M3.Ω, M3.ΩΛ/M3.Ω])
```
Notice that the variables vary "more nicely" with respect to $b$ than $t$, making the solver happier and avoiding numerical problems.
