using ModelingToolkit, OrdinaryDiffEq, Unitful
using Test
MT = ModelingToolkit
UMT = ModelingToolkit.UnitfulUnitCheck

@constants a = 1
@test_throws MT.ArgumentError @constants b

@independent_variables t
@variables x(t) w(t)
D = Differential(t)
eqs = [D(x) ~ a]
@named sys = ODESystem(eqs, t)
prob = ODEProblem(complete(sys), [0], [0.0, 1.0], [])
sol = solve(prob, Tsit5())

newsys = MT.eliminate_constants(sys)
@test isequal(equations(newsys), [D(x) ~ 1])

# Test structural_simplify substitutions & observed values
eqs = [D(x) ~ 1,
    w ~ a]
@named sys = ODESystem(eqs, t)
# Now eliminate the constants first
simp = structural_simplify(sys)
@test equations(simp) == [D(x) ~ 1.0]

#Constant with units
@constants β=1 [unit = u"m/s"]
UMT.get_unit(β)
@test MT.isconstant(β)
@independent_variables t [unit = u"s"]
@variables x(t) [unit = u"m"]
D = Differential(t)
eqs = [D(x) ~ β]
@named sys = ODESystem(eqs, t)
simp = structural_simplify(sys)

@test isempty(MT.collect_constants(nothing))

@testset "Issue#3044" begin
    @constants h = 1
    @parameters τ = 0.5 * h
    @variables x(MT.t_nounits) = h
    eqs = [MT.D_nounits(x) ~ (h - x) / τ]

    @mtkbuild fol_model = ODESystem(eqs, MT.t_nounits)

    prob = ODEProblem(fol_model, [], (0.0, 10.0))
    @test prob[x] ≈ 1
    @test prob.ps[τ] ≈ 0.5
end
