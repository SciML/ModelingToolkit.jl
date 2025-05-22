using ModelingToolkit, OrdinaryDiffEq, Unitful
using Test
MT = ModelingToolkit
UMT = ModelingToolkit.UnitfulUnitCheck

@constants a = 1
@test isconstant(a)
@test !istunable(a)

@independent_variables t
@variables x(t) w(t)
D = Differential(t)
eqs = [D(x) ~ a]
@named sys = System(eqs, t)
prob = ODEProblem(complete(sys), [0], [0.0, 1.0], [])
sol = solve(prob, Tsit5())

# Test structural_simplify substitutions & observed values
eqs = [D(x) ~ 1,
    w ~ a]
@named sys = System(eqs, t)
# Now eliminate the constants first
simp = structural_simplify(sys)
@test equations(simp) == [D(x) ~ 1.0]

#Constant with units
@constants β=1 [unit = u"m/s"]
UMT.get_unit(β)
@test MT.isconstant(β)
@test !MT.istunable(β)
@independent_variables t [unit = u"s"]
@variables x(t) [unit = u"m"]
D = Differential(t)
eqs = [D(x) ~ β]
@named sys = System(eqs, t)
simp = structural_simplify(sys)

@testset "Issue#3044" begin
    @constants h
    @parameters τ = 0.5 * h
    @variables x(MT.t_nounits) = h
    eqs = [MT.D_nounits(x) ~ (h - x) / τ]

    @mtkcompile fol_model = System(eqs, MT.t_nounits)

    prob = ODEProblem(fol_model, [], (0.0, 10.0), [h => 1])
    @test prob[x] ≈ 1
    @test prob.ps[τ] ≈ 0.5
end
