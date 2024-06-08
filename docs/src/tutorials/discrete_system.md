# (Experimental) Modeling Discrete Systems

In this example, we will use the new [`DiscreteSystem`](@ref) API
to create an SIR model.

```@example discrete
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using OrdinaryDiffEq: solve, FunctionMap

@inline function rate_to_proportion(r, t)
    1 - exp(-r * t)
end
@parameters c δt β γ
@constants h = 1
@variables S(t) I(t) R(t)
k = ShiftIndex(t)
infection = rate_to_proportion(
    β * c * I(k - 1) / (S(k - 1) * h + I(k - 1) + R(k - 1)), δt * h) * S(k - 1)
recovery = rate_to_proportion(γ * h, δt) * I(k - 1)

# Equations
eqs = [S(k) ~ S(k - 1) - infection * h,
    I(k) ~ I(k - 1) + infection - recovery,
    R(k) ~ R(k - 1) + recovery]
@mtkbuild sys = DiscreteSystem(eqs, t)

u0 = [S(k - 1) => 990.0, I(k - 1) => 10.0, R(k - 1) => 0.0]
p = [β => 0.05, c => 10.0, γ => 0.25, δt => 0.1]
tspan = (0.0, 100.0)
prob = DiscreteProblem(sys, u0, tspan, p)
sol = solve(prob, FunctionMap())
```

All shifts must be non-positive, i.e., discrete-time variables may only be indexed at index
`k, k-1, k-2, ...`. If default values are provided, they are treated as the value of the
variable at the previous timestep. For example, consider the following system to generate
the Fibonacci series:

```@example discrete
@variables x(t) = 1.0
@mtkbuild sys = DiscreteSystem([x ~ x(k - 1) + x(k - 2)], t)
```

The "default value" here should be interpreted as the value of `x` at all past timesteps.
For example, here `x(k-1)` and `x(k-2)` will be `1.0`, and the initial value of `x(k)` will
thus be `2.0`. During problem construction, the _past_ value of a variable should be
provided. For example, providing `[x => 1.0]` while constructing this problem will error.
Provide `[x(k-1) => 1.0]` instead. Note that values provided during problem construction
_do not_ apply to the entire history. Hence, if `[x(k-1) => 2.0]` is provided, the value of
`x(k-2)` will still be `1.0`.
