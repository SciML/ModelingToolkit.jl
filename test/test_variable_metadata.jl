using ModelingToolkit

# Bounds
@variables u [bounds = (-1, 1)]
@test getbounds(u) == (-1, 1)
@test hasbounds(u)
@test ModelingToolkit.dump_variable_metadata(u).bounds == (-1, 1)

@variables y
@test !hasbounds(y)
@test !haskey(ModelingToolkit.dump_variable_metadata(y), :bounds)

# Guess
@variables y [guess = 0]
@test getguess(y) === 0
@test hasguess(y) === true
@test ModelingToolkit.dump_variable_metadata(y).guess == 0

@variables y
@test hasguess(y) === false
@test !haskey(ModelingToolkit.dump_variable_metadata(y), :guess)

# Disturbance
@variables u [disturbance = true]
@test isdisturbance(u)
@test ModelingToolkit.dump_variable_metadata(u).disturbance

@variables y
@test !isdisturbance(y)
@test !haskey(ModelingToolkit.dump_variable_metadata(y), :disturbance)

# Tunable
@parameters u [tunable = true]
@test istunable(u)
@test ModelingToolkit.dump_variable_metadata(u).tunable

@parameters u2 [tunable = false]
@test !istunable(u2)
@test !ModelingToolkit.dump_variable_metadata(u2).tunable

@parameters y
@test istunable(y)
@test ModelingToolkit.dump_variable_metadata(y).tunable

# Distributions
struct FakeNormal end
d = FakeNormal()
@parameters u [dist = d]
@test hasdist(u)
@test getdist(u) == d
@test ModelingToolkit.dump_variable_metadata(u).dist == d

@parameters y
@test !hasdist(y)
@test !haskey(ModelingToolkit.dump_variable_metadata(y), :dist)

## System interface
@parameters t
Dₜ = Differential(t)
@variables x(t)=0 [bounds = (-10, 10)] u(t)=0 [input = true] y(t)=0 [output = true]
@parameters T [bounds = (0, Inf)]
@parameters k [tunable = true, bounds = (0, Inf)]
@parameters k2 [tunable = false]
eqs = [Dₜ(x) ~ (-k2 * x + k * u) / T
       y ~ x]
sys = ODESystem(eqs, t, name = :tunable_first_order)
unk_meta = ModelingToolkit.dump_unknowns(sys)
@test length(unk_meta) == 3
@test all(iszero, meta.default for meta in unk_meta)
param_meta = ModelingToolkit.dump_parameters(sys)
@test length(param_meta) == 3
@test all(!haskey(meta, :default) for meta in param_meta)

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

b = getbounds(sys, unknowns(sys))
@test b[x] == (-10, 10)

p = tunable_parameters(sys, default = false)
sp = Set(p)
@test k ∈ sp
@test T ∉ sp
@test k2 ∉ sp
@test length(p) == 1

## Descriptions
@variables u [description = "This is my input"]
@test getdescription(u) == "This is my input"
@test hasdescription(u)
@test ModelingToolkit.dump_variable_metadata(u).desc == "This is my input"

@variables u
@test getdescription(u) == ""
@test !hasdescription(u)
@test !haskey(ModelingToolkit.dump_variable_metadata(u), :desc)

@parameters t
@variables u(t) [description = "A short description of u"]
@parameters p [description = "A description of p"]
@named sys = ODESystem([u ~ p], t)

@test_nowarn show(stdout, "text/plain", sys)
