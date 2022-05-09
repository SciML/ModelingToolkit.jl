using ModelingToolkit

# Bounds
@variables u [bounds=(-1,1)]
@test getbounds(u) == (-1, 1)
@test hasbounds(u)

@variables y
@test !hasbounds(y)


# Disturbance
@variables u [disturbance=true]
@test isdisturbance(u)

@variables y
@test !isdisturbance(y)


# Tunable
@parameters u [tunable=true]
@test istunable(u)

@parameters y
@test !istunable(y)

# Distributions
struct FakeNormal end
d = FakeNormal()
@parameters u [dist=d]
@test hasdist(u)
@test getdist(u) == d

@parameters y
@test !hasdist(y)

## System interface
@parameters t
Dₜ = Differential(t)
@variables x(t)=0 u(t)=0 [input=true] y(t)=0 [output=true]
@parameters T [tunable = true, bounds = (0, Inf)]
@parameters k [tunable = true, bounds = (0, Inf)]
@parameters k2
eqs = [
    Dₜ(x) ~ (-k2*x + k*u) / T
    y ~ x
]
sys = ODESystem(eqs, t, name=:tunable_first_order)

p = tunable_parameters(sys)
sp = Set(p)
@test k ∈ sp
@test T ∈ sp
@test k2 ∉ sp
@test length(p) == 2

lb, ub = getbounds(p)
@test lb == [0,0]
@test ub == [Inf, Inf]

b = getbounds(sys)
@test b[T] == (0, Inf)

p = tunable_parameters(sys, default=true)
sp = Set(p)
@test k ∈ sp
@test T ∈ sp
@test k2 ∈ sp
@test length(p) == 3