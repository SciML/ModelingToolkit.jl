using ModelingToolkit

# Bounds
@variables u [bounds = (-1, 1)]
@test getbounds(u) == (-1, 1)
@test hasbounds(u)

@variables y
@test !hasbounds(y)

# Guess
@variables y [guess = 0]
@test getguess(y) === 0
@test hasguess(y) === true

@variables y
@test hasguess(y) === false

# Disturbance
@variables u [disturbance = true]
@test isdisturbance(u)

@variables y
@test !isdisturbance(y)

# Tunable
@parameters u [tunable = true]
@test istunable(u)

@parameters y
@test !istunable(y)

# Distributions
struct FakeNormal end
d = FakeNormal()
@parameters u [dist = d]
@test hasdist(u)
@test getdist(u) == d

@parameters y
@test !hasdist(y)

## System interface
@parameters t
Dₜ = Differential(t)
@variables x(t)=0 [bounds = (-10, 10)] u(t)=0 [input = true] y(t)=0 [output = true]
@parameters T [tunable = true, bounds = (0, Inf)]
@parameters k [tunable = true, bounds = (0, Inf)]
@parameters k2
eqs = [Dₜ(x) ~ (-k2 * x + k * u) / T
    y ~ x]
sys = ODESystem(eqs, t, name = :tunable_first_order)

p = tunable_parameters(sys)
sp = Set(p)
@test k ∈ sp
@test T ∈ sp
@test k2 ∉ sp
@test length(p) == 2

lb, ub = getbounds(p)
@test lb == [0, 0]
@test ub == [Inf, Inf]

b = getbounds(sys)
@test b[T] == (0, Inf)

b = getbounds(sys, states(sys))
@test b[x] == (-10, 10)

p = tunable_parameters(sys, default = true)
sp = Set(p)
@test k ∈ sp
@test T ∈ sp
@test k2 ∈ sp
@test length(p) == 3

## Descriptions
@variables u [description = "This is my input"]
@test getdescription(u) == "This is my input"
@test hasdescription(u)

@variables u
@test getdescription(u) == ""
@test !hasdescription(u)

@parameters t
@variables u(t) [description = "A short description of u"]
@parameters p [description = "A description of p"]
@named sys = ODESystem([u ~ p], t)

@test_nowarn show(stdout, "text/plain", sys)

@testset "binary" begin
    @parameters t
    @variables u(t) [binary = true]
    @parameters p [binary = true]
    @test isbinaryvar(u)
    @test isbinaryvar(p)
end

@testset "integer" begin
    @parameters t
    @variables u(t) [integer = true]
    @parameters p [integer = true]
    @test isintegervar(u)
    @test isintegervar(p)
end
