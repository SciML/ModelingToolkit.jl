using ModelingToolkit, StaticArrays, LinearAlgebra
using DiffEqBase
using Test

# Define some variables
@independent_variables t
@parameters σ ρ β
@variables x(t) y(t) z(t)
D = Differential(t)

eqs = [D(x) ~ σ * (y - x),
    D(y) ~ -z - y,
    D(z) ~ y - β * z]

@test ModelingToolkit.islinear(@named sys = ODESystem(eqs, t))

eqs2 = [D(x) ~ σ * (y - x),
    D(y) ~ -z - 1 / y,
    D(z) ~ y - β * z]

@test !ModelingToolkit.islinear(@named sys = ODESystem(eqs2, t))

eqs3 = [D(x) ~ σ * (y - x),
    D(y) ~ -z - y,
    D(z) ~ y - β * z + 1]

@test ModelingToolkit.isaffine(@named sys = ODESystem(eqs, t))
