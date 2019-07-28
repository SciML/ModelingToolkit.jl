using ModelingToolkit, StaticArrays, LinearAlgebra
using DiffEqBase
using Test

# Define some variables
@parameters t σ ρ β
@variables x y z

eqs = [σ*(y-x),
       x*(ρ-z)-y,
       x*y - β*z]

∂ = ModelingToolkit.jacobian(eqs,[x,y,z])
for i in 1:3
    ∇ = ModelingToolkit.gradient(eqs[i],[x,y,z])
    @test isequal(∂[i,:],∇)
end
@test all(isequal.(ModelingToolkit.gradient(eqs[1],[x,y,z]),[σ * -1,σ,0]))
@test all(isequal.(ModelingToolkit.hessian(eqs[1],[x,y,z]),0))
