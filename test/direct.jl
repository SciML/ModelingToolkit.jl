using ModelingToolkit, StaticArrays, LinearAlgebra
using DiffEqBase
using Test

# Calculus
@parameters t σ ρ β
@variables x y z

eqs = [σ*(y-x),
       x*(ρ-z)-y,
       x*y - β*z]

simpexpr = [
   :(σ * (y - x))
   :(x * (ρ - z) - y)
   :(x * y - β * z)
   ]

for i in 1:3
   @test ModelingToolkit.simplified_expr.(eqs)[i] == simpexpr[i]
   @test ModelingToolkit.simplified_expr.(eqs)[i] == simpexpr[i]
end

∂ = ModelingToolkit.jacobian(eqs,[x,y,z])
for i in 1:3
    ∇ = ModelingToolkit.gradient(eqs[i],[x,y,z])
    @test isequal(∂[i,:],∇)
end
@test all(isequal.(ModelingToolkit.gradient(eqs[1],[x,y,z]),[σ * -1,σ,0]))
@test all(isequal.(ModelingToolkit.hessian(eqs[1],[x,y,z]),0))

# Function building

@parameters σ() ρ() β()
@variables x y z
eqs = [σ*(y-x),
       x*(ρ-z)-y,
       x*y - β*z]
f = eval(ModelingToolkit.build_function(eqs,[x,y,z],[σ,ρ,β]))
out = [1.0,2,3]
o1 = f([1.0,2,3],[1.0,2,3])
f(out,[1.0,2,3],[1.0,2,3])
@test all(o1 .== out)

function test_worldage()
   @parameters σ() ρ() β()
   @variables x y z
   eqs = [σ*(y-x),
          x*(ρ-z)-y,
          x*y - β*z]
   f, f_iip = ModelingToolkit.build_function(eqs,[x,y,z],[σ,ρ,β],(),ModelingToolkit.simplified_expr,Val{false})
   out = [1.0,2,3]
   o1 = f([1.0,2,3],[1.0,2,3])
   f_iip(out,[1.0,2,3],[1.0,2,3])
end
test_worldage()

mac = @I begin
   @parameters σ() ρ() β()
   @variables x() y() z()

   eqs = [σ*(y-x),
          x*(ρ-z)-y,
          x*y - β*z]
   ModelingToolkit.build_function(eqs,[x,y,z],[σ,ρ,β])
end
f = @ICompile
out = [1.0,2,3]
o1 = f([1.0,2,3],[1.0,2,3])
f(out,[1.0,2,3],[1.0,2,3])
@test all(o1 .== out)

mac = @I begin
   @parameters σ ρ β
   @variables x y z

   eqs = [σ*(y-x),
          x*(ρ-z)-y,
          x*y - β*z]
   ∂ = ModelingToolkit.jacobian(eqs,[x,y,z])
   ModelingToolkit.build_function(∂,[x,y,z],[σ,ρ,β])
end
f = @ICompile
out = zeros(3,3)
o1 = f([1.0,2,3],[1.0,2,3])
f(out,[1.0,2,3],[1.0,2,3])
@test all(out .== o1)

## No parameters
@variables x y z
eqs = [(y-x)^2,
       x*(x-z)-y,
       x*y - y*z]
f = eval(ModelingToolkit.build_function(eqs,[x,y,z]))
out = zeros(3)
o1 = f([1.0,2,3])
f(out,[1.0,2,3])
@test all(out .== o1)

function test_worldage()
   @variables x y z
   eqs = [(y-x)^2,
          x*(x-z)-y,
          x*y - y*z]
   f, f_iip = ModelingToolkit.build_function(eqs,[x,y,z],(),(),ModelingToolkit.simplified_expr,Val{false})
   out = zeros(3)
   o1 = f([1.0,2,3])
   f_iip(out,[1.0,2,3])
end
test_worldage()

mac = @I begin
   @variables x y z
   eqs = [(y-x)^2,
          x*(x-z)-y,
          x*y - y*z]
   ModelingToolkit.build_function(eqs,[x,y,z])
end
f = @ICompile
out = zeros(3)
o1 = f([1.0,2,3])
f(out,[1.0,2,3])
@test all(out .== o1)
