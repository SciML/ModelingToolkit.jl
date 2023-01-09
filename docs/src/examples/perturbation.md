# [Symbolic-Numeric Perturbation Theory for ODEs](@id perturb_diff)

## Prelims

In the previous tutorial, [Mixed Symbolic-Numeric Perturbation Theory](https://symbolics.juliasymbolics.org/stable/examples/perturbation), we discussed how to solve algebraic equations using **Symbolics.jl**. Here, our goal is to extend the method to differential equations. First, we import the following helper functions that were introduced in [Mixed Symbolic/Numerical Methods for Perturbation Theory - Algebraic Equations](@ref perturb_alg):

```julia
using Symbolics, SymbolicUtils

def_taylor(x, ps) = sum([a*x^i for (i,a) in enumerate(ps)])
def_taylor(x, ps, p₀) = p₀ + def_taylor(x, ps)

function collect_powers(eq, x, ns; max_power=100)
    eq = substitute(expand(eq), Dict(x^j => 0 for j=last(ns)+1:max_power))

    eqs = []
    for i in ns
        powers = Dict(x^j => (i==j ? 1 : 0) for j=1:last(ns))
        push!(eqs, substitute(eq, powers))
    end
    eqs
end

function solve_coef(eqs, ps)
    vals = Dict()

    for i = 1:length(ps)
        eq = substitute(eqs[i], vals)
        vals[ps[i]] = Symbolics.solve_for(eq ~ 0, ps[i])
    end
    vals
end
```

## The Trajectory of a Ball!

In the first two examples, we applied the perturbation method to algebraic problems. However, the main power of the perturbation method is to solve differential equations (usually ODEs, but also occasionally PDEs). Surprisingly, the main procedure developed to solve algebraic problems works well for differential equations. In fact, we will use the same two helper functions, `collect_powers` and `solve_coef`. The main difference is in the way we expand the dependent variables. For algebraic problems, the coefficients of $\epsilon$ are constants; whereas, for differential equations, they are functions of the dependent variable (usually time).

As the first ODE example, we have chosen a simple and well-behaved problem, which is a variation of a standard first-year physics problem: what is the trajectory of an object (say, a ball, or a rocket) thrown vertically at velocity $v$ from the surface of a planet? Assuming a constant acceleration of gravity, $g$, every burgeoning physicist knows the answer: $x(t) = x(0) + vt - \frac{1}{2}gt^2$. However, what happens if $g$ is not constant? Specifically, $g$ is inversely proportional to the distant from the center of the planet. If $v$ is large and the projectile travels a large fraction of the radius of the planet, the assumption of constant gravity does not hold anymore. However, unless $v$ is large compared to the escape velocity, the correction is usually small. After simplifications and change of variables to dimensionless, the problem becomes

```math
  \ddot{x}(t) = -\frac{1}{(1 + \epsilon x(t))^2}
```

with the initial conditions $x(0) = 0$, and $\dot{x}(0) = 1$. Note that for $\epsilon = 0$, this equation transforms back to the standard one. Let's start with defining the variables

```julia
n = 3
@variables ϵ t y[1:n](t) ∂∂y[1:n](t)
```

Next, we define $x$.

```julia
x = def_taylor(ϵ, y[3:end], y[2])
```

We need the second derivative of `x`. It may seem that we can do this using `Differential(t)`; however, this operation needs to wait for a few steps because we need to manipulate the differentials as separate variables. Instead, we define dummy variables `∂∂y` as the placeholder for the second derivatives and define

```julia
∂∂x = def_taylor(ϵ, ∂∂y[3:end], ∂∂y[2])
```

as the second derivative of `x`. After rearrangement, our governing equation is $\ddot{x}(t)(1 + \epsilon x(t))^{-2} + 1 = 0$, or

```julia
eq = ∂∂x * (1 + ϵ*x)^2 + 1
```

The next two steps are the same as the ones for algebraic equations (note that we pass `1:n` to `collect_powers` because the zeroth order term is needed here)

```julia
eqs = collect_powers(eq, ϵ, 1:n)
```

and,

```julia
vals = solve_coef(eqs, ∂∂y)
```

Our system of ODEs is forming. Now is the time to convert `∂∂`s to the correct **Symbolics.jl** form by substitution:

```julia
D = Differential(t)
subs = Dict(∂∂y[i] => D(D(y[i])) for i in eachindex(y))
eqs = [substitute(first(v), subs) ~ substitute(last(v), subs) for v in vals]
```

We are nearly there! From this point on, the rest is standard ODE solving procedures. Potentially, we can use a symbolic ODE solver to find a closed form solution to this problem. However, **Symbolics.jl** currently does not support this functionality. Instead, we solve the problem numerically. We form an `ODESystem`, lower the order (convert second derivatives to first), generate an `ODEProblem` (after passing the correct initial conditions), and, finally, solve it.

```julia
using ModelingToolkit, DifferentialEquations

@named sys = ODESystem(eqs, t)
sys = structural_simplify(sys)
states(sys)
```

```julia
# the initial conditions
# everything is zero except the initial velocity
u0 = zeros(2n+2)
u0[3] = 1.0   # y₀ˍt

prob = ODEProblem(sys, u0, (0, 3.0))
sol = solve(prob; dtmax=0.01)
```

Finally, we calculate the solution to the problem as a function of `ϵ` by substituting the solution to the ODE system back into the defining equation for `x`. Note that `𝜀` is a number, compared to `ϵ`, which is a symbolic variable.

```julia
X = 𝜀 -> sum([𝜀^(i-1) * sol[y[i]] for i in eachindex(y)])
```

Using `X`, we can plot the trajectory for a range of $𝜀$s.

```julia
using Plots

plot(sol.t, hcat([X(𝜀) for 𝜀 = 0.0:0.1:0.5]...))
```

As expected, as `𝜀` becomes larger (meaning the gravity is less with altitude), the object goes higher and stays up for a longer duration. Of course, we could have solved the problem directly using as ODE solver. One of the benefits of the perturbation method is that we need to run the ODE solver only once and then can just calculate the answer for different values of `𝜀`; whereas, if we had used the direct method, we would need to run the solver once for each value of `𝜀`.

## A Weakly Nonlinear Oscillator

For the next example, we have chosen a simple example from a very important class of problems, the nonlinear oscillators. As we will see, perturbation theory has difficulty providing a good solution to this problem, but the process is instructive. This example closely follows the chapter 7.6 of *Nonlinear Dynamics and Chaos* by Steven Strogatz.

The goal is to solve $\ddot{x} + 2\epsilon\dot{x} + x = 0$, where the dot signifies time-derivatives and the initial conditions are $x(0) = 0$ and $\dot{x}(0) = 1$. If $\epsilon = 0$, the problem reduces to the simple linear harmonic oscillator with the exact solution $x(t) = \sin(t)$. We follow the same steps as the previous example.

```julia
n = 3
@variables ϵ t y[1:n](t) ∂y[1:n] ∂∂y[1:n]
x = def_taylor(ϵ, y[3:end], y[2])
∂x = def_taylor(ϵ, ∂y[3:end], ∂y[2])
∂∂x = def_taylor(ϵ, ∂∂y[3:end], ∂∂y[2])
```

This time we also need the first derivative terms. Continuing,

```julia
eq = ∂∂x + 2*ϵ*∂x + x
eqs = collect_powers(eq, ϵ, 0:n)
vals = solve_coef(eqs, ∂∂y)
```

Next, we need to replace `∂`s and `∂∂`s with their **Symbolics.jl** counterparts:

```julia
D = Differential(t)
subs1 = Dict(∂y[i] => D(y[i]) for i in eachindex(y))
subs2 = Dict(∂∂y[i] => D(D(y[i])) for i in eachindex(y))
subs = subs1 ∪ subs2
eqs = [substitute(first(v), subs) ~ substitute(last(v), subs) for v in vals]
```

We continue with converting 'eqs' to an `ODEProblem`, solving it, and finally plot the results against the exact solution to the original problem, which is $x(t, \epsilon) = (1 - \epsilon)^{-1/2} e^{-\epsilon t} \sin((1- \epsilon^2)^{1/2}t)$,

```julia
@named sys = ODESystem(eqs, t)
sys = structural_simplify(sys)
```

```julia
# the initial conditions
u0 = zeros(2n+2)
u0[3] = 1.0   # y₀ˍt

prob = ODEProblem(sys, u0, (0, 50.0))
sol = solve(prob; dtmax=0.01)

X = 𝜀 -> sum([𝜀^(i-1) * sol[y[i]] for i in eachindex(y)])
T = sol.t
Y = 𝜀 -> exp.(-𝜀*T) .* sin.(sqrt(1 - 𝜀^2)*T) / sqrt(1 - 𝜀^2)    # exact solution

plot(sol.t, [Y(0.1), X(0.1)])
```

The figure is similar to Figure 7.6.2 in *Nonlinear Dynamics and Chaos*. The two curves fit well for the first couple of cycles, but then the perturbation method curve diverges from the true solution. The main reason is that the problem has two or more time-scales that introduce secular terms in the solution. One solution is to explicitly account for the two time scales and use an analytic method called *two-timing*.
