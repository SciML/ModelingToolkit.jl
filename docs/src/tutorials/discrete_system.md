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
infection = rate_to_proportion(β * c * I(k-1) / (S(k-1) * h + I(k-1) + R(k-1)), δt * h) * S(k-1)
recovery = rate_to_proportion(γ * h, δt) * I(k-1)

# Equations
eqs = [S(k) ~ S(k-1) - infection * h,
    I(k) ~ I(k-1) + infection - recovery,
    R(k) ~ R(k-1) + recovery]
@mtkbuild sys = DiscreteSystem(eqs, t)

u0 = [S(k - 1) => 990.0, I(k - 1) => 10.0, R(k - 1) => 0.0]
p = [β => 0.05, c => 10.0, γ => 0.25, δt => 0.1]
tspan = (0.0, 100.0)
prob = DiscreteProblem(sys, u0, tspan, p)
sol = solve(prob, FunctionMap())
```

All shifts must be negative. If default values are provided, they are treated as the value
for the variable at the previous timestep. For example, consider the following system to
generate the Fibonacci series:

```@example discrete
@variables x(t) = 1.0
@mtkbuild sys = DiscreteSystem([x ~ x(k-1) + x(k-2)], t)
```

Note that the default value is treated as the initial value of `x(k-1)`. The value for
`x(k-2)` must be provided during problem construction.
