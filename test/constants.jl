using ModelingToolkit, OrdinaryDiffEq, Unitful
using Test
MT = ModelingToolkit

@constants a = 1
@test_throws MT.ArgumentError @constants b

@variables t x(t) w(t)
D = Differential(t)
eqs = [D(x) ~ a]
@named sys = ODESystem(eqs)
prob = ODEProblem(sys, [0], [0.0, 1.0], [])
sol = solve(prob, Tsit5())

# Test structural_simplify substitutions & observed values
eqs = [D(x) ~ 1,
    w ~ a]
@named sys = ODESystem(eqs)
simp = structural_simplify(sys);
@test isequal(simp.substitutions.subs[1], eqs[2])
@test isequal(equations(simp)[1], eqs[1])
prob = ODEProblem(simp, [0], [0.0, 1.0], [])
sol = solve(prob, Tsit5())
@test sol[w][1] == 1

#Constant with units
@constants β=1 [unit = u"m/s"]
MT.get_unit(β)
@test MT.isconstant(β)
@variables t [unit = u"s"] x(t) [unit = u"m"]
D = Differential(t)
eqs = [D(x) ~ β]
sys = ODESystem(eqs, name = :sys)
