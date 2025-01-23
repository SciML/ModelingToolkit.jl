using ModelingToolkit
using DynamicQuantities

# Bounds
@variables u [bounds = (-1, 1)]
@test getbounds(u) == (-1, 1)
@test hasbounds(u)
@test ModelingToolkit.dump_variable_metadata(u).bounds == (-1, 1)

@variables y
@test !hasbounds(y)
@test !haskey(ModelingToolkit.dump_variable_metadata(y), :bounds)

@variables y[1:3]
@test !hasbounds(y)
@test getbounds(y)[1] == [-Inf, -Inf, -Inf]
@test getbounds(y)[2] == [Inf, Inf, Inf]
for i in eachindex(y)
    @test !hasbounds(y[i])
    b = getbounds(y[i])
    @test b[1] == -Inf && b[2] == Inf
end

@variables y[1:3] [bounds = (-1, 1)]
@test hasbounds(y)
@test getbounds(y)[1] == -ones(3)
@test getbounds(y)[2] == ones(3)
for i in eachindex(y)
    @test hasbounds(y[i])
    b = getbounds(y[i])
    @test b[1] == -1.0 && b[2] == 1.0
end
@test getbounds(y[1:2])[1] == -ones(2)
@test getbounds(y[1:2])[2] == ones(2)

@variables y[1:2, 1:2] [bounds = (-1, [1.0 Inf; 2.0 3.0])]
@test hasbounds(y)
@test getbounds(y)[1] == [-1 -1; -1 -1]
@test getbounds(y)[2] == [1.0 Inf; 2.0 3.0]
for i in eachindex(y)
    @test hasbounds(y[i])
    b = getbounds(y[i])
    @test b[1] == -1 && b[2] == [1.0 Inf; 2.0 3.0][i]
end

@variables y[1:2] [bounds = (-Inf, [1.0, Inf])]
@test hasbounds(y)
@test getbounds(y)[1] == [-Inf, -Inf]
@test getbounds(y)[2] == [1.0, Inf]
@test hasbounds(y[1])
@test getbounds(y[1]) == (-Inf, 1.0)
@test !hasbounds(y[2])
@test getbounds(y[2]) == (-Inf, Inf)

# Guess
@variables y [guess = 0]
@test getguess(y) === 0
@test hasguess(y) === true
@test ModelingToolkit.dump_variable_metadata(y).guess == 0

# Default
@variables y = 0
@test ModelingToolkit.getdefault(y) === 0
@test ModelingToolkit.hasdefault(y) === true
@test ModelingToolkit.dump_variable_metadata(y).default == 0

# Issue#2653
@variables y[1:3] [guess = ones(3)]
@test getguess(y) == ones(3)
@test hasguess(y) === true
@test ModelingToolkit.dump_variable_metadata(y).guess == ones(3)

for i in 1:3
    @test getguess(y[i]) == 1.0
    @test hasguess(y[i]) === true
    @test ModelingToolkit.dump_variable_metadata(y[i]).guess == 1.0
end

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
@independent_variables t
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

@independent_variables t
@variables u(t) [description = "A short description of u"]
@parameters p [description = "A description of p"]
@named sys = ODESystem([u ~ p], t)

@test_nowarn show(stdout, "text/plain", sys)

# Defaults, guesses overridden by system, parameter dependencies
@variables x(t)=1.0 y(t) [guess = 1.0]
@parameters p=2.0 q
@named sys = ODESystem(Equation[], t, [x, y], [p]; defaults = Dict(x => 2.0, p => 3.0),
    guesses = Dict(y => 2.0), parameter_dependencies = [q => 2p])
unks_meta = ModelingToolkit.dump_unknowns(sys)
unks_meta = Dict([ModelingToolkit.getname(meta.var) => meta for meta in unks_meta])
@test unks_meta[:x].default == 2.0
@test unks_meta[:y].guess == 2.0
params_meta = ModelingToolkit.dump_parameters(sys)
params_meta = Dict([ModelingToolkit.getname(meta.var) => meta for meta in params_meta])
@test params_meta[:p].default == 3.0
@test isequal(params_meta[:q].dependency, 2p)

# Connect
@variables x [connect = Flow]
@test hasconnect(x)
@test getconnect(x) == Flow
@test ModelingToolkit.dump_variable_metadata(x).connect == Flow
x = ModelingToolkit.setconnect(x, ModelingToolkit.Stream)
@test getconnect(x) == ModelingToolkit.Stream

struct BadConnect end
@test_throws Exception ModelingToolkit.setconnect(x, BadConnect)

# Unit
@variables x [unit = u"s"]
@test hasunit(x)
@test getunit(x) == u"s"
@test ModelingToolkit.dump_variable_metadata(x).unit == u"s"

# Misc data
@variables x [misc = [:good]]
@test hasmisc(x)
@test getmisc(x) == [:good]
x = ModelingToolkit.setmisc(x, "okay")
@test getmisc(x) == "okay"

# Variable Type
@variables x
@test ModelingToolkit.getvariabletype(x) == ModelingToolkit.VARIABLE
@test ModelingToolkit.dump_variable_metadata(x).variable_type == ModelingToolkit.VARIABLE
@test ModelingToolkit.dump_variable_metadata(x).variable_source == :variables
x = ModelingToolkit.toparam(x)
@test ModelingToolkit.getvariabletype(x) == ModelingToolkit.PARAMETER
@test ModelingToolkit.dump_variable_metadata(x).variable_source == :variables

@parameters y
@test ModelingToolkit.getvariabletype(y) == ModelingToolkit.PARAMETER

@brownian z
@test ModelingToolkit.getvariabletype(z) == ModelingToolkit.BROWNIAN
