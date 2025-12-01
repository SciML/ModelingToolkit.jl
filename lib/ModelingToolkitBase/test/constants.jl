using ModelingToolkitBase, OrdinaryDiffEq, DynamicQuantities
using Test
MT = ModelingToolkitBase

@constants a = 1
@test isconstant(a)
@test !istunable(a)

@independent_variables t
@variables x(t) w(t)
D = Differential(t)
eqs = [D(x) ~ a]
@named sys = System(eqs, t)
prob = ODEProblem(complete(sys), [x => 0], [0.0, 1.0])
sol = solve(prob, Tsit5())

# Test mtkcompile substitutions & observed values
eqs = [D(x) ~ 1,
    w ~ a]
@named sys = System(eqs, t)
# Now eliminate the constants first
simp = mtkcompile(sys)
if @isdefined(ModelingToolkit)
    @test equations(simp) == [D(x) ~ 1.0]
else
    @test equations(simp) == [D(x) ~ 1.0, 0 ~ a-w]
end

#Constant with units
@constants β=1 [unit = u"m/s"]
MT.get_unit(β)
@test MT.isconstant(β)
@test !MT.istunable(β)
@independent_variables t [unit = u"s"]
@variables x(t) [unit = u"m"]
D = Differential(t)
eqs = [D(x) ~ β]
@named sys = System(eqs, t)
simp = mtkcompile(sys)

@testset "Issue#3044" begin
    @constants h
    @parameters τ = 0.5 * h
    @variables x(MT.t_nounits) = h
    eqs = [MT.D_nounits(x) ~ (h - x) / τ]

    @mtkcompile fol_model = System(eqs, MT.t_nounits)

    prob = ODEProblem(fol_model, [h => 1], (0.0, 10.0))
    @test prob[x] ≈ 1
    @test prob.ps[τ] ≈ 0.5
end
