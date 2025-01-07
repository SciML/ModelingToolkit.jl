using ModelingToolkit, FMI, FMIZoo, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkit as MTK

const FMU_DIR = joinpath(@__DIR__, "fmus")

@testset "Standalone pendulum model" begin
    fmu = loadFMU("SpringPendulum1D", "Dymola", "2022x"; type = :ME)
    truesol = FMI.simulate(
        fmu, (0.0, 8.0); saveat = 0.0:0.1:8.0, recordValues = ["mass.s", "mass.v"])

    function test_no_inputs_outputs(sys)
        for var in unknowns(sys)
            @test !MTK.isinput(var)
            @test !MTK.isoutput(var)
        end
    end
    @testset "v2, ME" begin
        fmu = loadFMU("SpringPendulum1D", "Dymola", "2022x"; type = :ME)
        @mtkbuild sys = MTK.FMIComponent(Val(2); fmu, type = :ME)
        test_no_inputs_outputs(sys)
        prob = ODEProblem{true, SciMLBase.FullSpecialize}(
            sys, [sys.mass__s => 0.5, sys.mass__v => 0.0], (0.0, 8.0))
        sol = solve(prob, Tsit5(); reltol = 1e-8, abstol = 1e-8)
        @test SciMLBase.successful_retcode(sol)

        @test sol(0.0:0.1:8.0;
            idxs = [sys.mass__s, sys.mass__v]).u≈collect.(truesol.values.saveval) atol=1e-4
        # repeated solve works
        @test_nowarn solve(prob, Tsit5())
    end
    @testset "v2, CS" begin
        fmu = loadFMU("SpringPendulum1D", "Dymola", "2022x"; type = :CS)
        @named inner = MTK.FMIComponent(
            Val(2); fmu, communication_step_size = 0.001, type = :CS)
        @variables x(t) = 1.0
        @mtkbuild sys = ODESystem([D(x) ~ x], t; systems = [inner])
        test_no_inputs_outputs(sys)

        prob = ODEProblem{true, SciMLBase.FullSpecialize}(
            sys, [sys.inner.mass__s => 0.5, sys.inner.mass__v => 0.0], (0.0, 8.0))
        sol = solve(prob, Tsit5(); reltol = 1e-8, abstol = 1e-8)
        @test SciMLBase.successful_retcode(sol)

        @test sol(0.0:0.1:8.0;
            idxs = [sys.inner.mass__s, sys.inner.mass__v]).u≈collect.(truesol.values.saveval) rtol=1e-2
    end

    fmu = loadFMU("SpringPendulum1D", "Dymola", "2023x", "3.0"; type = :ME)
    truesol = FMI.simulate(
        fmu, (0.0, 8.0); saveat = 0.0:0.1:8.0, recordValues = ["mass.s", "mass.v"])
    @testset "v3, ME" begin
        fmu = loadFMU("SpringPendulum1D", "Dymola", "2023x", "3.0"; type = :ME)
        @mtkbuild sys = MTK.FMIComponent(Val(3); fmu, type = :ME)
        test_no_inputs_outputs(sys)
        prob = ODEProblem{true, SciMLBase.FullSpecialize}(
            sys, [sys.mass__s => 0.5, sys.mass__v => 0.0], (0.0, 8.0))
        sol = solve(prob, Tsit5(); reltol = 1e-8, abstol = 1e-8)
        @test SciMLBase.successful_retcode(sol)

        @test sol(0.0:0.1:8.0;
            idxs = [sys.mass__s, sys.mass__v]).u≈collect.(truesol.values.saveval) atol=1e-4
        # repeated solve works
        @test_nowarn solve(prob, Tsit5())
    end
    @testset "v3, CS" begin
        fmu = loadFMU("SpringPendulum1D", "Dymola", "2023x", "3.0"; type = :CS)
        @named inner = MTK.FMIComponent(
            Val(3); fmu, communication_step_size = 0.001, type = :CS)
        @variables x(t) = 1.0
        @mtkbuild sys = ODESystem([D(x) ~ x], t; systems = [inner])
        test_no_inputs_outputs(sys)

        prob = ODEProblem{true, SciMLBase.FullSpecialize}(
            sys, [sys.inner.mass__s => 0.5, sys.inner.mass__v => 0.0], (0.0, 8.0))
        sol = solve(prob, Tsit5(); reltol = 1e-8, abstol = 1e-8)
        @test SciMLBase.successful_retcode(sol)

        @test sol(0.0:0.1:8.0;
            idxs = [sys.inner.mass__s, sys.inner.mass__v]).u≈collect.(truesol.values.saveval) rtol=1e-2
    end
end

@testset "IO Model" begin
    @testset "v2, ME" begin
        fmu = loadFMU(joinpath(FMU_DIR, "SimpleAdder.fmu"); type = :ME)
        @named adder = MTK.FMIComponent(Val(2); fmu, type = :ME)
        @test MTK.isinput(adder.a)
        @test MTK.isinput(adder.b)
        @test MTK.isoutput(adder.out)
        @test !MTK.isinput(adder.c) && !MTK.isoutput(adder.c)

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
    @testset "v2, CS" begin
        fmu = loadFMU(joinpath(FMU_DIR, "SimpleAdder.fmu"); type = :CS)
        @named adder = MTK.FMIComponent(
            Val(2); fmu, type = :CS, communication_step_size = 0.001)
        @test MTK.isinput(adder.a)
        @test MTK.isinput(adder.b)
        @test MTK.isoutput(adder.out)
        @test !MTK.isinput(adder.c) && !MTK.isoutput(adder.c)
        @variables a(t) b(t) c(t) [guess = 1.0]
        @mtkbuild sys = ODESystem(
            [adder.a ~ a, adder.b ~ b, D(a) ~ t,
                D(b) ~ adder.out + adder.c, c^2 ~ adder.out + adder.value],
            t;
            systems = [adder])

        # c will be solved for by initialization
        # this tests that initialization also works with FMUs
        prob = ODEProblem(sys, [sys.adder.c => 1.0, sys.a => 1.0, sys.b => 1.0],
            (0.0, 1.0); use_scc = false)
        sol = solve(prob, Rodas5P(autodiff = false))
        @test SciMLBase.successful_retcode(sol)
    end

    @testset "v3, ME" begin
        fmu = loadFMU(joinpath(FMU_DIR, "StateSpace.fmu"); type = :ME)
        @named sspace = MTK.FMIComponent(Val(3); fmu, type = :ME)
        @test MTK.isinput(sspace.u)
        @test MTK.isoutput(sspace.y)
        @test !MTK.isinput(sspace.x) && !MTK.isoutput(sspace.x)
        @variables u(t)=1.0 x(t)=1.0 y(t) [guess = 1.0]
        @mtkbuild sys = ODESystem(
            [sspace.u ~ u, D(u) ~ t, D(x) ~ sspace.x + sspace.y, y^2 ~ sspace.y + x], t;
            systems = [sspace]
        )

        prob = ODEProblem(sys, [sys.sspace.x => 1.0], (0.0, 1.0); use_scc = false)
        sol = solve(prob, Rodas5P(autodiff = false))
        @test SciMLBase.successful_retcode(sol)
    end

    @testset "v3, CS" begin
        fmu = loadFMU(joinpath(FMU_DIR, "StateSpace.fmu"); type = :CS)
        @named sspace = MTK.FMIComponent(
            Val(3); fmu, communication_step_size = 1e-3, type = :CS)
        @test MTK.isinput(sspace.u)
        @test MTK.isoutput(sspace.y)
        @test !MTK.isinput(sspace.x) && !MTK.isoutput(sspace.x)
        @variables u(t)=1.0 x(t)=1.0 y(t) [guess = 1.0]
        @mtkbuild sys = ODESystem(
            [sspace.u ~ u, D(u) ~ t, D(x) ~ sspace.x + sspace.y, y^2 ~ sspace.y + x], t;
            systems = [sspace]
        )

        prob = ODEProblem(sys, [sys.sspace.x => 1.0], (0.0, 1.0); use_scc = false)
        sol = solve(prob, Rodas5P(autodiff = false))
        @test SciMLBase.successful_retcode(sol)
    end
end
