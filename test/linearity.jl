using ModelingToolkit, StaticArrays, LinearAlgebra
using DiffEqBase
using Test

# Define some variables
@parameters t σ ρ β
@variables x(t) y(t) z(t)
@derivatives D'~t

# Define a differential equation
eqs = [D(x) ~ σ*(y-x),
       D(y) ~ -z-y,
       D(z) ~ y - β*z]

sys = ODESystem(eqs)

@test ModelingToolkit.islinear(sys)
