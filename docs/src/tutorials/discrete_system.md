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
infection = rate_to_proportion(β * c * I / (S * h + I + R), δt * h) * S
recovery = rate_to_proportion(γ * h, δt) * I

# Equations
eqs = [S(k + 1) ~ S - infection * h,
    I(k + 1) ~ I + infection - recovery,
    R(k + 1) ~ R + recovery]
@mtkbuild sys = DiscreteSystem(eqs, t)

u0 = [S => 990.0, I => 10.0, R => 0.0]
p = [β => 0.05, c => 10.0, γ => 0.25, δt => 0.1]
tspan = (0.0, 100.0)
prob = DiscreteProblem(sys, u0, tspan, p)
sol = solve(prob, FunctionMap())
```
