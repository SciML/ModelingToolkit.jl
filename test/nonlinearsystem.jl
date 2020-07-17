using ModelingToolkit, StaticArrays, LinearAlgebra
using DiffEqBase, SparseArrays
using Test

canonequal(a, b) = isequal(simplify(a), simplify(b))

# Define some variables
@parameters t σ ρ β
@variables x y z

function test_nlsys_inference(name, sys, vs, ps)
    @testset "NonlinearSystem construction: $name" begin
        @test Set(states(sys))     == Set(convert.(Variable,vs))
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

# Now nonlinear system with only variables
@variables x y z
@parameters σ ρ β

# Define a nonlinear system
eqs = [0 ~ σ*(y-x),
       0 ~ x*(ρ-z)-y,
       0 ~ x*y - β*z]
ns = NonlinearSystem(eqs, [x,y,z], [σ,ρ,β])
jac = calculate_jacobian(ns)
@testset "nlsys jacobian" begin
    @test canonequal(jac[1,1], σ * -1)
    @test canonequal(jac[1,2], σ)
    @test canonequal(jac[1,3], 0)
    @test canonequal(jac[2,1], ρ - z)
    @test canonequal(jac[2,2], -1)
    @test canonequal(jac[2,3], x * -1)
    @test canonequal(jac[3,1], y)
    @test canonequal(jac[3,2], x)
    @test canonequal(jac[3,3], -1 * β)
end
nlsys_func = generate_function(ns, [x,y,z], [σ,ρ,β])
jac_func = generate_jacobian(ns)
f = @eval eval(nlsys_func)

# Intermediate calculations
a = y - x
# Define a nonlinear system
eqs = [0 ~ σ*a,
       0 ~ x*(ρ-z)-y,
       0 ~ x*y - β*z]
ns = NonlinearSystem(eqs, [x,y,z], [σ,ρ,β])
nlsys_func = generate_function(ns, [x,y,z], [σ,ρ,β])
jac = calculate_jacobian(ns)

@test ModelingToolkit.jacobian_sparsity(ns).colptr == sparse(jac).colptr
@test ModelingToolkit.jacobian_sparsity(ns).rowval == sparse(jac).rowval

jac = generate_jacobian(ns)
