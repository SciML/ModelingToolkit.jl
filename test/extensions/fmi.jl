using ModelingToolkit, FMI, FMIZoo, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkit as MTK

@testset "Standalone pendulum model" begin
    @testset "v2, ME" begin
        fmu = loadFMU("SpringPendulum1D", "Dymola", "2022x"; type = :ME)
        @mtkbuild sys = MTK.FMIComponent(Val(2), Val(:ME); fmu)
        prob = ODEProblem{true, SciMLBase.FullSpecialize}(
            sys, [sys.mass__s => 0.5, sys.mass__v => 0.0], (0.0, 8.0))
        sol = solve(prob, Tsit5(); reltol = 1e-8, abstol = 1e-8)
        @test SciMLBase.successful_retcode(sol)

        truesol = FMI.simulate(fmu, (0.0, 8.0); saveat = 0.0:0.1:8.0)
        @test sol(0.0:0.1:8.0).u≈truesol.states.u atol=1e-4
        # repeated solve works
        @test_nowarn solve(prob, Tsit5())
    end
end

@testset "IO Model" begin
    fmu = loadFMU("./fmus/SimpleAdder.fmu"; type = :ME)
    @named adder = MTK.FMIComponent(Val(2), Val(:ME); fmu)
    @variables a(t) b(t) c(t) [guess = 1.0]
    @mtkbuild sys = ODESystem(
        [adder.a ~ a, adder.b ~ b, D(a) ~ t,
            D(b) ~ adder.out + adder.c, c^2 ~ adder.out + adder.value],
        t;
        systems = [adder])

    # c will be solved for by initialization
    # this tests that initialization also works with FMUs
    prob = ODEProblem(sys, [sys.adder.c => 1.0, sys.a => 1.0, sys.b => 1.0], (0.0, 1.0))
    sol = solve(prob, Rodas5P(autodiff = false))
    @test SciMLBase.successful_retcode(sol)
end
