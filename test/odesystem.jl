using ModelingToolkit, StaticArrays, LinearAlgebra
using OrdinaryDiffEq
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

ModelingToolkit.simplified_expr.(eqs)[1]
:(derivative(x(t), t) = σ * (y(t) - x(t))).args
de = ODESystem(eqs)

generate_function(de)

function _clean(O::Operation)
    @assert isa(O.op, Variable)
    return O.op
end
_clean(x::Variable) = x
function test_diffeq_inference(name, sys, iv, dvs, ps)
    @testset "ODESystem construction: $name" begin
        @test independent_variable(sys) ==  _clean(iv)
        @test Set(states(sys))          == Set(_clean.(dvs))
        @test Set(parameters(sys))      == Set(_clean.(ps))
    end
end

test_diffeq_inference("standard", de, t, [x, y, z], [ρ, σ, β])
generate_function(de, [x,y,z], [σ,ρ,β])
jac_expr = generate_jacobian(de)
jac = calculate_jacobian(de)
jacfun = eval(jac_expr[2])
# iip
f = ODEFunction(de, [x,y,z], [σ,ρ,β])
Wfact, Wfact_t = ModelingToolkit.calculate_factorized_W(de)
fw, fwt = map(eval, ModelingToolkit.generate_factorized_W(de))
du = zeros(3)
u  = collect(1:3)
p  = collect(4:6)
f(du, u, p, 0.1)
@test du == [4, 0, -16]
J = zeros(3, 3)
jacfun(J, u, p, t)
FW = zeros(3, 3)
FWt = zeros(3, 3)
eval(fw[2])(FW, u, p, 0.2, 0.1)
eval(fwt[2])(FWt, u, p, 0.2, 0.1)
# oop
f = ODEFunction(de, [x,y,z], [σ,ρ,β])
fw, fwt = map(eval, ModelingToolkit.generate_factorized_W(de))
du = @SArray zeros(3)
u  = SVector(1:3...)
p  = SVector(4:6...)
@test f(u, p, 0.1) === @SArray [4, 0, -16]
Sfw = eval(fw[1])(u, p, 0.2, 0.1)
@test Sfw.L ≈ UnitLowerTriangular(FW)
@test Sfw.U ≈ UpperTriangular(FW)
sol = Sfw \ @SArray ones(3)
@test sol isa SArray
@test sol ≈ -(I - 0.2*J)\ones(3)
Sfw_t = eval(fwt[1])(u, p, 0.2, 0.1)
@test Sfw_t.L ≈ UnitLowerTriangular(FWt)
@test Sfw_t.U ≈ UpperTriangular(FWt)
sol = Sfw_t \ @SArray ones(3)
@test sol isa SArray
@test sol ≈ -(I/0.2 - J)\ones(3)

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y*t,
       D(z) ~ x*y - β*z]
de = ODESystem(eqs)
ModelingToolkit.calculate_tgrad(de)

tgrad_oop, tgrad_iip = eval.(ModelingToolkit.generate_tgrad(de))
@test tgrad_oop(u,p,t) == [0.0,-u[2],0.0]
du = zeros(3)
tgrad_iip(du,u,p,t)
@test du == [0.0,-u[2],0.0]

@testset "time-varying parameters" begin
    @parameters σ′(t-1)
    eqs = [D(x) ~ σ′*(y-x),
           D(y) ~ x*(ρ-z)-y,
           D(z) ~ x*y - β*z]
    de = ODESystem(eqs)
    test_diffeq_inference("global iv-varying", de, t, (x, y, z), (σ′, ρ, β))
    @test begin
        f = eval(generate_function(de, [x,y,z], [σ′,ρ,β])[2])
        du = [0.0,0.0,0.0]
        f(du, [1.0,2.0,3.0], [x->x+7,2,3], 5.0)
        du ≈ [11, -3, -7]
    end

    @parameters σ(..)
    eqs = [D(x) ~ σ(t-1)*(y-x),
           D(y) ~ x*(ρ-z)-y,
           D(z) ~ x*y - β*z]
    de = ODESystem(eqs)
    test_diffeq_inference("single internal iv-varying", de, t, (x, y, z), (σ, ρ, β))
    @test begin
        f = eval(generate_function(de, [x,y,z], [σ,ρ,β])[2])
        du = [0.0,0.0,0.0]
        f(du, [1.0,2.0,3.0], [x->x+7,2,3], 5.0)
        du ≈ [11, -3, -7]
    end

    eqs = [D(x) ~ x + 10σ(t-1) + 100σ(t-2) + 1000σ(t^2)]
    de = ODESystem(eqs)
    test_diffeq_inference("many internal iv-varying", de, t, (x,), (σ,))
    @test begin
        f = eval(generate_function(de, [x], [σ])[2])
        du = [0.0]
        f(du, [1.0], [t -> t + 2], 5.0)
        du ≈ [27561]
    end
end

# Conversion to first-order ODEs #17
@derivatives D3'''~t
@derivatives D2''~t
@variables u(t) uˍtt(t) uˍt(t) xˍt(t)
eqs = [D3(u) ~ 2(D2(u)) + D(u) + D(x) + 1
       D2(x) ~ D(x) + 2]
de = ODESystem(eqs)
de1 = ode_order_lowering(de)
lowered_eqs = [D(uˍtt) ~ 2uˍtt + uˍt + xˍt + 1
               D(xˍt)  ~ xˍt + 2
               D(uˍt)  ~ uˍtt
               D(u)    ~ uˍt
               D(x)    ~ xˍt]

@test de1 == ODESystem(lowered_eqs)

# issue #219
@test de1.states == [ModelingToolkit.var_from_nested_derivative(eq.lhs)[1] for eq in de1.eqs] == ODESystem(lowered_eqs).states

test_diffeq_inference("first-order transform", de1, t, [uˍtt, xˍt, uˍt, u, x], [])
du = zeros(5)
ODEFunction(de1, [uˍtt, xˍt, uˍt, u, x], [])(du, ones(5), nothing, 0.1)
@test du == [5.0, 3.0, 1.0, 1.0, 1.0]

# Internal calculations
a = y - x
eqs = [D(x) ~ σ*a,
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]
de = ODESystem(eqs)
generate_function(de, [x,y,z], [σ,ρ,β])
jac = calculate_jacobian(de)
f = ODEFunction(de, [x,y,z], [σ,ρ,β])

@derivatives D'~t
@parameters A B C
_x = y / C
eqs = [D(x) ~ -A*x,
       D(y) ~ A*x - B*_x]
de = ODESystem(eqs)
@test begin
    f = eval(generate_function(de, [x,y], [A,B,C])[2])
    du = [0.0,0.0]
    f(du, [1.0,2.0], [1,2,3], 0.0)
    du ≈ [-1, -1/3]
    f = eval(generate_function(de, [x,y], [A,B,C])[1])
    du ≈ f([1.0,2.0], [1,2,3], 0.0)
end

function lotka(u,p,t)
  x = u[1]
  y = u[2]
  [p[1]*x - p[2]*x*y,
  -p[3]*y + p[4]*x*y]
end

prob = ODEProblem(ODEFunction{false}(lotka),[1.0,1.0],(0.0,1.0),[1.5,1.0,3.0,1.0])
de = modelingtoolkitize(prob)
ODEFunction(de)(similar(prob.u0), prob.u0, prob.p, 0.1)

function lotka(du,u,p,t)
  x = u[1]
  y = u[2]
  du[1] = p[1]*x - p[2]*x*y
  du[2] = -p[3]*y + p[4]*x*y
end

prob = ODEProblem(lotka,[1.0,1.0],(0.0,1.0),[1.5,1.0,3.0,1.0])

de = modelingtoolkitize(prob)
ODEFunction(de)(similar(prob.u0), prob.u0, prob.p, 0.1)

# automatic state detection for DAEs
@parameters t k₁ k₂ k₃
@variables y₁(t) y₂(t) y₃(t)
@derivatives D'~t
# reorder the system just to be a little spicier
eqs = [D(y₁) ~ -k₁*y₁+k₃*y₂*y₃,
       0     ~  y₁ + y₂ + y₃ - 1,
       D(y₂) ~  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃]
sys = ODESystem(eqs)
u0 = [y₁ => 1.0,
      y₂ => 0.0,
      y₃ => 0.0]
p  = [k₁ => 0.04,
      k₂ => 3e7,
      k₃ => 1e4]
tspan = (0.0,100000.0)
prob1 = ODEProblem(sys,u0,tspan,p)
prob2 = ODEProblem(sys,u0,tspan,p,jac=true)
prob3 = ODEProblem(sys,u0,tspan,p,Wfact=true,Wfact_t=true)
for (prob, atol) in [(prob1, 1e-12), (prob2, 1e-12), (prob3, 1e-12)]
    sol = solve(prob, Rodas5())
    @test all(x->≈(sum(x), 1.0, atol=atol), sol.u)
end
