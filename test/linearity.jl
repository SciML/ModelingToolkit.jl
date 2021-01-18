using ModelingToolkit, StaticArrays, LinearAlgebra
using DiffEqBase
using Test

# Define some variables
@parameters t σ ρ β
@variables x(t) y(t) z(t)
D = Differential(t)

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ -z-y,
       D(z) ~ y - β*z]

@test ModelingToolkit.islinear(ODESystem(eqs))

eqs2 = [D(x) ~ σ*(y-x),
       D(y) ~ -z-1/y,
       D(z) ~ y - β*z]

@test !ModelingToolkit.islinear(ODESystem(eqs2))
