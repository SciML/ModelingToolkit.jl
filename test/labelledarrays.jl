using ModelingToolkit, StaticArrays, LinearAlgebra, LabelledArrays
using DiffEqBase
using Test

# Define some variables
@parameters t σ ρ β
@variables x(t) y(t) z(t)
@derivatives D'~t

# Define a differential equation
eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]

de = ODESystem(eqs)
f = ODEFunction(de, [x,y,z], [σ,ρ,β])

a = @SVector [1.0,2.0,3.0]
b = SLVector(x=1.0,y=2.0,z=3.0)
c = [1.0,2.0,3.0]
p = SLVector(σ=10.0,ρ=26.0,β=8/3)
@test f(a,p,0.0) isa SVector
@test typeof(f(b,p,0.0)) <: SLArray
@test f(c,p,0.0) isa Vector
@test f(a,p,0.0) == f(b,p,0.0)
@test f(a,p,0.0) == f(c,p,0.0)
