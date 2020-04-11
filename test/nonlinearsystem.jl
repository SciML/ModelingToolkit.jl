using ModelingToolkit, StaticArrays, LinearAlgebra
using DiffEqBase
using Test

# Define some variables
@parameters t σ ρ β
@variables x(t) y(t) z(t)

function test_nlsys_inference(name, sys, vs, ps)
    @testset "NonlinearSystem construction: $name" begin
        @test Set(states(sys))     == Set(vs)
        @test Set(parameters(sys)) == Set(convert.(Variable,ps))
    end
end

# Define a nonlinear system
eqs = [0 ~ σ*(y-x),
       0 ~ x*(ρ-z)-y,
       0 ~ x*y - β*z]
ns = NonlinearSystem(eqs, [x,y,z], [σ,ρ,β])
test_nlsys_inference("standard", ns, (x, y, z), (σ, ρ, β))
@test begin
    f = eval(generate_function(ns, [x,y,z], [σ,ρ,β])[2])
    du = [0.0, 0.0, 0.0]
    f(du, [1,2,3], [1,2,3])
    du ≈ [1, -3, -7]
end
