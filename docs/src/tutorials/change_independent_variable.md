# Changing the independent variable of ODEs

Ordinary differential equations describe the rate of change of some dependent variables with respect to one independent variable.
For the modeler it is often most natural to write down the equations with a particular independent variable, say time $t$.
However, in many cases there are good reasons for changing the independent variable:

 1. One may want $y(x)$ as a function of $x$ instead of $(x(t), y(t))$ as a function of $t$

 2. Some differential equations vary more nicely (e.g. less stiff) with respect to one independent variable than another.
 3. It can reduce the number of equations that must be solved (e.g. $y(x)$ is one equation, while $(x(t), y(t))$ are two).

To manually change the independent variable of an ODE, one must rewrite all equations in terms of a new variable and transform differentials with the chain rule.
This is mechanical and error-prone.
ModelingToolkit provides the utility function [`change_independent_variable`](@ref) that automates this process.

## 1. Get one dependent variable as a function of another

Consider a projectile shot with some initial velocity in a vertical gravitational field with constant horizontal velocity.

```@example changeivar
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
@variables x(t) y(t)
@parameters g=9.81 v # gravitational acceleration and horizontal velocity
eqs = [D(D(y)) ~ -g, D(x) ~ v]
initialization_eqs = [D(x) ~ D(y)] # 45° initial angle
M1 = System(eqs, t; initialization_eqs, name = :M)
M1s = mtkcompile(M1)
@assert length(equations(M1s)) == 3 # hide
M1s # hide
```

This is the standard parametrization that arises naturally from kinematics and Newton's laws.
It expresses the position $(x(t), y(t))$ as a function of time $t$.
But suppose we want to determine whether the projectile hits a target 10 meters away.
There are at least three ways of answering this:

  - Solve the ODE for $(x(t), y(t))$ and use a callback to terminate when $x$ reaches 10 meters, and evaluate $y$ at the final time.
  - Solve the ODE for $(x(t), y(t))$ and use root finding to find the time when $x$ reaches 10 meters, and evaluate $y$ at that time.
  - Solve the ODE for $y(x)$ and evaluate it at 10 meters.

We will demonstrate the last method by changing the independent variable from $t$ to $x$.
This transformation is well-defined for any non-zero horizontal velocity $v$, so $x$ and $t$ are one-to-one.

```@example changeivar
M2 = change_independent_variable(M1, x)
M2s = mtkcompile(M2; allow_symbolic = true)
# a sanity test on the 10 x/y variables that are accessible to the user # hide
@assert allequal([x, M1s.x]) # hide
@assert allequal([M2.x, M2s.x]) # hide
@assert allequal([y, M1s.y]) # hide
@assert allunique([M1.x, M1.y, M2.y, M2s.y]) # hide
@assert length(equations(M2s)) == 2 # hide
M2s # display this # hide
```

The derivatives are now with respect to the new independent variable $x$, which can be accessed with `M2.x`.

!!! info
    
    At this point `x`, `M1.x`, `M1s.x`, `M2.x`, `M2s.x` are *three* different variables.
    Meanwhile `y`, `M1.y`, `M1s.y`, `M2.y` and `M2s.y` are *four* different variables.
    It can be instructive to inspect these yourself to see their subtle differences.

Notice how the number of equations has also decreased from three to two, as $\mathrm{d}x/\mathrm{d}t$ has been turned into an observed equation.
It is straightforward to evolve the ODE for 10 meters and plot the resulting trajectory $y(x)$:

```@example changeivar
using OrdinaryDiffEq, Plots
prob = ODEProblem(M2s, [M2s.y => 0.0, v => 8.0], [0.0, 10.0]) # throw 10 meters
sol = solve(prob, Tsit5())
plot(sol; idxs = M2.y) # must index by M2.y = y(x); not M1.y = y(t)!
```

!!! tip "Usage tips"
    
    Look up the documentation of [`change_independent_variable`](@ref) for tips on how to use it.
    
    For example, if you also need $t(x)$, you can tell it to add a differential equation for the old independent variable in terms of the new one using the [inverse function rule](https://en.wikipedia.org/wiki/Inverse_function_rule) (here $\mathrm{d}t/\mathrm{d}x = 1 / (\mathrm{d}x/\mathrm{d}t)$). If you know an analytical expression between the independent variables (here $t = x/v$), you can also pass it directly to the function to avoid the extra differential equation.

## 2. Alleviating stiffness by changing to logarithmic time

In cosmology, the [Friedmann equations](https://en.wikipedia.org/wiki/Friedmann_equations) describe the expansion of the universe.
In terms of conformal time $t$, they can be written

```@example changeivar
@variables a(t) Ω(t)
a = GlobalScope(a) # global var needed by all species
function species(w; kw...)
    eqs = [D(Ω) ~ -3(1 + w) * D(a) / a * Ω]
    return System(eqs, t, [Ω], []; kw...)
end
@named r = species(1 // 3) # radiation
@named m = species(0) # matter
@named Λ = species(-1) # dark energy / cosmological constant
eqs = [Ω ~ r.Ω + m.Ω + Λ.Ω, D(a) ~ √(Ω) * a^2]
initialization_eqs = [Λ.Ω + r.Ω + m.Ω ~ 1]
M1 = System(eqs, t, [Ω, a], []; initialization_eqs, name = :M)
M1 = compose(M1, r, m, Λ)
M1s = mtkcompile(M1)
```

Of course, we can solve this ODE as it is:

```@example changeivar
prob = ODEProblem(M1s, [M1s.a => 1.0, M1s.r.Ω => 5e-5, M1s.m.Ω => 0.3], (0.0, -3.3), [])
sol = solve(prob, Tsit5(); reltol = 1e-5)
@assert Symbol(sol.retcode) == :Unstable # surrounding text assumes this was unstable # hide
plot(sol, idxs = [M1.a, M1.r.Ω / M1.Ω, M1.m.Ω / M1.Ω, M1.Λ.Ω / M1.Ω])
```

But the solver becomes unstable due to stiffness.
Also notice the interesting dynamics taking place towards the end of the integration (in the early universe), which gets compressed into a very small time interval.
These ODEs would benefit from being defined with respect to a logarithmic "time" that better captures the evolution of the universe through *orders of magnitude* of time, rather than linear time.

It is therefore common to write these ODEs in terms of $b = \ln a$.
To do this, we will change the independent variable in two stages; first from $t$ to $a$, and then from $a$ to $b$.
Notice that $\mathrm{d}a/\mathrm{d}t > 0$ provided that $\Omega > 0$, and $\mathrm{d}b/\mathrm{d}a > 0$, so the transformation is well-defined since $t \leftrightarrow a \leftrightarrow b$ are one-to-one.
First, we transform from $t$ to $a$:

```@example changeivar
M2 = change_independent_variable(M1, M1.a)
@assert !ModelingToolkit.isautonomous(M2) # hide
M2 # hide
```

Unlike the original, notice that this system is *non-autonomous* because the independent variable $a$ appears explicitly in the equations!
This means that to change the independent variable from $a$ to $b$, we must provide not only the rate of change relation $db(a)/da = \exp(-b)$, but *also* the equation $a(b) = \exp(b)$ so $a$ can be eliminated in favor of $b$:

```@example changeivar
a = M2.a # get independent variable of M2
Da = Differential(a)
@variables b(a)
M3 = change_independent_variable(M2, b, [Da(b) ~ exp(-b), a ~ exp(b)])
```

We can now solve and plot the ODE in terms of $b$:

```@example changeivar
M3s = mtkcompile(M3; allow_symbolic = true)
prob = ODEProblem(M3s, [M3s.r.Ω => 5e-5, M3s.m.Ω => 0.3], (0, -15), [])
sol = solve(prob, Tsit5())
@assert Symbol(sol.retcode) == :Success # surrounding text assumes the solution was successful # hide
plot(sol, idxs = [M3.r.Ω / M3.Ω, M3.m.Ω / M3.Ω, M3.Λ.Ω / M3.Ω])
```

Notice that the variables vary "more nicely" with respect to $b$ than $t$, making the solver happier and avoiding numerical problems.
