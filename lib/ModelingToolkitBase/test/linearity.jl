using ModelingToolkitBase, StaticArrays, LinearAlgebra
using DiffEqBase
using Test

# Define some variables
@independent_variables t
@parameters σ ρ β
@variables x(t) y(t) z(t)
D = Differential(t)

eqs = [
    D(x) ~ σ * (y - x),
    D(y) ~ -z - y,
    D(z) ~ y - β * z,
]

@test ModelingToolkitBase.islinear(@named sys = System(eqs, t))

eqs2 = [
    D(x) ~ σ * (y - x),
    D(y) ~ -z - 1 / y,
    D(z) ~ y - β * z,
]

@test !ModelingToolkitBase.islinear(@named sys = System(eqs2, t))

eqs3 = [
    D(x) ~ σ * (y - x),
    D(y) ~ -z - y,
    D(z) ~ y - β * z + 1,
]

@test ModelingToolkitBase.isaffine(@named sys = System(eqs, t))
