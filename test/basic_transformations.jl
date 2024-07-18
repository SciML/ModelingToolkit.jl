using ModelingToolkit, OrdinaryDiffEq, Test

@independent_variables t
@parameters α β γ δ
@variables x(t) y(t)
D = Differential(t)

eqs = [D(x) ~ α * x - β * x * y,
    D(y) ~ -δ * y + γ * x * y]

sys = ODESystem(eqs)

u0 = [x => 1.0,
    y => 1.0]

p = [α => 1.5,
    β => 1.0,
    δ => 3.0,
    γ => 1.0]

tspan = (0.0, 10.0)
prob = ODEProblem(sys, u0, tspan, p)
sol = solve(prob, Tsit5())

sys2 = liouville_transform(sys)
@variables trJ

u0 = [x => 1.0,
    y => 1.0,
    trJ => 1.0]

prob = ODEProblem(sys2, u0, tspan, p, jac = true)
sol = solve(prob, Tsit5())
@test sol[end, end] ≈ 1.0742818931017244
