using ModelingToolkit, OrdinaryDiffEq, Test, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters σ ρ β
@variables x(t) y(t) z(t) k(t)

eqs = [D(D(x)) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z]

@named sys′ = ODESystem(eqs, t)
sys = ode_order_lowering(sys′)

eqs2 = [0 ~ x * y - k,
    D(D(x)) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z]
@named sys2 = ODESystem(eqs2, t, [x, y, z, k], parameters(sys′))
sys2 = ode_order_lowering(sys2)
# test equation/variable ordering
ModelingToolkit.calculate_massmatrix(sys2) == Diagonal([1, 1, 1, 1, 0])

u0 = [D(x) => 2.0,
    x => 1.0,
    y => 0.0,
    z => 0.0]

p = [σ => 28.0,
    ρ => 10.0,
    β => 8 / 3]

tspan = (0.0, 100.0)

sys = complete(sys)
prob = ODEProblem(sys, u0, tspan, p, jac = true)
probexpr = ODEProblemExpr(sys, u0, tspan, p, jac = true)
sol = solve(prob, Tsit5())
solexpr = solve(eval(prob), Tsit5())
@test all(x -> x == 0, Array(sol - solexpr))
#using Plots; plot(sol,idxs=(:x,:y))

@parameters σ ρ β
@variables x(t) y(t) z(t)

eqs = [D(x) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z]

lorenz1 = ODESystem(eqs, t, name = :lorenz1)
lorenz2 = ODESystem(eqs, t, name = :lorenz2)

@variables α(t)
@parameters γ
connections = [0 ~ lorenz1.x + lorenz2.y + α * γ]
@named connected = ODESystem(connections, t, [α], [γ], systems = [lorenz1, lorenz2])
connected = complete(connected)
u0 = [lorenz1.x => 1.0,
    lorenz1.y => 0.0,
    lorenz1.z => 0.0,
    lorenz2.x => 0.0,
    lorenz2.y => 1.0,
    lorenz2.z => 0.0,
    α => 2.0]

p = [lorenz1.σ => 10.0,
    lorenz1.ρ => 28.0,
    lorenz1.β => 8 / 3,
    lorenz2.σ => 10.0,
    lorenz2.ρ => 28.0,
    lorenz2.β => 8 / 3,
    γ => 2.0]

tspan = (0.0, 100.0)
prob = ODEProblem(connected, u0, tspan, p)
sol = solve(prob, Rodas5())
@test maximum(sol[2, :] + sol[6, :] + 2sol[1, :]) < 1e-12
#using Plots; plot(sol,idxs=(:α,Symbol(lorenz1.x),Symbol(lorenz2.y)))
