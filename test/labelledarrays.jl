using ModelingToolkit, StaticArrays, LinearAlgebra, LabelledArrays
using DiffEqBase, ForwardDiff
using Test

# Define some variables
@parameters t σ ρ β
@variables x(t) y(t) z(t)
D = Differential(t)

# Define a differential equation
eqs = [D(x) ~ σ*(y-x),
       D(y) ~ t*x*(ρ-z)-y,
       D(z) ~ x*y - β*z]

de = ODESystem(eqs)
ff = ODEFunction(de, [x,y,z], [σ,ρ,β], jac=true)

a = @SVector [1.0,2.0,3.0]
b = SLVector(x=1.0,y=2.0,z=3.0)
c = [1.0,2.0,3.0]
p = SLVector(σ=10.0,ρ=26.0,β=8/3)
@test ff(a,p,0.0) isa SVector
@test typeof(ff(b,p,0.0)) <: SLArray
@test ff(c,p,0.0) isa Vector
@test ff(a,p,0.0) == ff(b,p,0.0)
@test ff(a,p,0.0) == ff(c,p,0.0)

@test ff.jac(a,p,0.0) isa SMatrix
@test typeof(ff.jac(b,p,0.0)) <: SMatrix
@test ff.jac(c,p,0.0) isa Matrix
@test ff.jac(a,p,0.0) == ff.jac(b,p,0.0)
@test ff.jac(a,p,0.0) == ff.jac(c,p,0.0)

# Test similar_type
@test ff(b,p,ForwardDiff.Dual(0.0,1.0)) isa SLArray
d = LVector(x=1.0,y=2.0,z=3.0)
@test ff(d,p,ForwardDiff.Dual(0.0,1.0)) isa LArray
@test ff.jac(b,p,ForwardDiff.Dual(0.0,1.0)) isa SArray
@test eltype(ff.jac(b,p,ForwardDiff.Dual(0.0,1.0))) <: ForwardDiff.Dual
@test ff.jac(d,p,ForwardDiff.Dual(0.0,1.0)) isa Array
@inferred ff.jac(d,p,ForwardDiff.Dual(0.0,1.0))
@test eltype(ff.jac(d,p,ForwardDiff.Dual(0.0,1.0))) <: ForwardDiff.Dual
