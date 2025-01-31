using ModelingToolkit, FMI, FMIZoo, OrdinaryDiffEq, NonlinearSolve, SciMLBase
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

@mtkmodel SimpleAdder begin
    @variables begin
        a(t)
        b(t)
        c(t)
        out(t)
        out2(t)
    end
    @parameters begin
        value = 1.0
    end
    @equations begin
        out ~ a + b + value
        D(c) ~ out
        out2 ~ 2c
    end
end

@mtkmodel StateSpace begin
    @variables begin
        x(t)
        y(t)
        u(t)
    end
    @parameters begin
        A = 1.0
        B = 1.0
        C = 1.0
        _D = 1.0
    end
    @equations begin
        D(x) ~ A * x + B * u
        y ~ C * x + _D * u
    end
end

@testset "IO Model" begin
    function build_simple_adder(adder)
        @variables a(t) b(t) c(t) [guess = 1.0]
        @mtkbuild sys = ODESystem(
            [adder.a ~ a, adder.b ~ b, D(a) ~ t,
                D(b) ~ adder.out + adder.c, c^2 ~ adder.out + adder.value],
            t;
            systems = [adder])
        # c will be solved for by initialization
        # this tests that initialization also works with FMUs
        prob = ODEProblem(sys, [sys.adder.c => 2.0, sys.a => 1.0, sys.b => 1.0],
            (0.0, 1.0), [sys.adder.value => 2.0])
        return sys, prob
    end

    @named adder = SimpleAdder()
    truesys, trueprob = build_simple_adder(adder)
    truesol = solve(trueprob, abstol = 1e-8, reltol = 1e-8)
    @test SciMLBase.successful_retcode(truesol)

    @testset "v2, ME" begin
        fmu = loadFMU(joinpath(FMU_DIR, "SimpleAdder.fmu"); type = :ME)
        @named adder = MTK.FMIComponent(Val(2); fmu, type = :ME)
        @test MTK.isinput(adder.a)
        @test MTK.isinput(adder.b)
        @test MTK.isoutput(adder.out)
        @test MTK.isoutput(adder.out2)
        @test !MTK.isinput(adder.c) && !MTK.isoutput(adder.c)

        sys, prob = build_simple_adder(adder)
        sol = solve(prob, Rodas5P(autodiff = false), abstol = 1e-8, reltol = 1e-8)
        @test SciMLBase.successful_retcode(sol)

        @test truesol(sol.t;
            idxs = [truesys.a, truesys.b, truesys.c, truesys.adder.c]).u≈sol[[
            sys.a, sys.b, sys.c, sys.adder.c]] rtol=1e-7
    end
    @testset "v2, CS" begin
        fmu = loadFMU(joinpath(FMU_DIR, "SimpleAdder.fmu"); type = :CS)
        @named adder = MTK.FMIComponent(
            Val(2); fmu, type = :CS, communication_step_size = 1e-3,
            reinitializealg = BrownFullBasicInit())
        @test MTK.isinput(adder.a)
        @test MTK.isinput(adder.b)
        @test MTK.isoutput(adder.out)
        @test MTK.isoutput(adder.out2)
        @test !MTK.isinput(adder.c) && !MTK.isoutput(adder.c)

        sys, prob = build_simple_adder(adder)
        sol = solve(prob, Rodas5P(autodiff = false), abstol = 1e-8, reltol = 1e-8)
        @test SciMLBase.successful_retcode(sol)

        @test truesol(sol.t; idxs = [truesys.a, truesys.b, truesys.c]).u≈sol.u rtol=1e-2
        # sys.adder.c is a discrete variable
        @test truesol(sol.t; idxs = truesys.adder.c).u≈sol(sol.t; idxs = sys.adder.c).u rtol=1e-3
    end

    function build_sspace_model(sspace)
        @variables u(t)=1.0 x(t)=1.0 y(t) [guess = 1.0]
        @mtkbuild sys = ODESystem(
            [sspace.u ~ u, D(u) ~ t, D(x) ~ sspace.x + sspace.y, y^2 ~ sspace.y + sspace.x], t;
            systems = [sspace]
        )

        prob = ODEProblem(
            sys, [sys.sspace.x => 1.0], (0.0, 1.0), [sys.sspace.A => 2.0]; use_scc = false)
        return sys, prob
    end

    @named sspace = StateSpace()
    truesys, trueprob = build_sspace_model(sspace)
    truesol = solve(trueprob, abstol = 1e-8, reltol = 1e-8)
    @test SciMLBase.successful_retcode(truesol)

    @testset "v3, ME" begin
        fmu = loadFMU(joinpath(FMU_DIR, "StateSpace.fmu"); type = :ME)
        @named sspace = MTK.FMIComponent(Val(3); fmu, type = :ME)
        @test MTK.isinput(sspace.u)
        @test MTK.isoutput(sspace.y)
        @test !MTK.isinput(sspace.x) && !MTK.isoutput(sspace.x)

        sys, prob = build_sspace_model(sspace)
        sol = solve(prob, Rodas5P(autodiff = false); abstol = 1e-8, reltol = 1e-8)
        @test SciMLBase.successful_retcode(sol)

        @test truesol(sol.t;
            idxs = [truesys.u, truesys.x, truesys.y, truesys.sspace.x]).u≈sol[[
            sys.u, sys.x, sys.y, sys.sspace.x]] rtol=1e-7
    end

    @testset "v3, CS" begin
        fmu = loadFMU(joinpath(FMU_DIR, "StateSpace.fmu"); type = :CS)
        @named sspace = MTK.FMIComponent(
            Val(3); fmu, communication_step_size = 1e-4, type = :CS,
            reinitializealg = BrownFullBasicInit())
        @test MTK.isinput(sspace.u)
        @test MTK.isoutput(sspace.y)
        @test !MTK.isinput(sspace.x) && !MTK.isoutput(sspace.x)

        sys, prob = build_sspace_model(sspace)
        sol = solve(prob, Rodas5P(autodiff = false); abstol = 1e-8, reltol = 1e-8)
        @test SciMLBase.successful_retcode(sol)

        @test truesol(
            sol.t; idxs = [truesys.u, truesys.x, truesys.y]).u≈sol[[sys.u, sys.x, sys.y]] rtol=1e-2
        @test truesol(sol.t; idxs = truesys.sspace.x).u≈sol(sol.t; idxs = sys.sspace.x).u rtol=1e-2
    end
end

@testset "FMUs in a loop" begin
    function build_looped_adders(adder1, adder2)
        @variables x(t) = 1
        @mtkbuild sys = ODESystem(
            [D(x) ~ x, adder1.a ~ adder2.out2,
                adder2.a ~ adder1.out2, adder1.b ~ 1.0, adder2.b ~ 2.0],
            t;
            systems = [adder1, adder2])
        prob = ODEProblem(
            sys, [adder1.c => 1.0, adder2.c => 1.0, adder1.a => 2.0], (0.0, 1.0))
        return sys, prob
    end
    @named adder1 = SimpleAdder()
    @named adder2 = SimpleAdder()
    truesys, trueprob = build_looped_adders(adder1, adder2)
    truesol = solve(trueprob, Tsit5(), reltol = 1e-8)
    @test SciMLBase.successful_retcode(truesol)

    @testset "v2, ME" begin
        fmu = loadFMU(joinpath(FMU_DIR, "SimpleAdder.fmu"); type = :ME)
        @named adder1 = MTK.FMIComponent(Val(2); fmu, type = :ME)
        @named adder2 = MTK.FMIComponent(Val(2); fmu, type = :ME)
        sys, prob = build_looped_adders(adder1, adder2)
        sol = solve(prob, Rodas5P(autodiff = false); reltol = 1e-8)
        @test SciMLBase.successful_retcode(sol)
        @test truesol(sol.t;
            idxs = [truesys.adder1.c, truesys.adder2.c]).u≈sol(
            sol.t; idxs = [sys.adder1.c, sys.adder2.c]).u rtol=1e-7
    end
    @testset "v2, CS" begin
        fmu = loadFMU(joinpath(FMU_DIR, "SimpleAdder.fmu"); type = :CS)
        @named adder1 = MTK.FMIComponent(
            Val(2); fmu, type = :CS, communication_step_size = 1e-3)
        @named adder2 = MTK.FMIComponent(
            Val(2); fmu, type = :CS, communication_step_size = 1e-3)
        sys, prob = build_looped_adders(adder1, adder2)
        sol = solve(prob,
            Tsit5();
            reltol = 1e-8,
            initializealg = SciMLBase.OverrideInit(nlsolve = FastShortcutNLLSPolyalg(autodiff = AutoFiniteDiff())))
        @test truesol(sol.t;
            idxs = [truesys.adder1.c, truesys.adder2.c]).u≈sol(
            sol.t; idxs = [sys.adder1.c, sys.adder2.c]).u rtol=1e-3
    end

    function build_looped_sspace(sspace1, sspace2)
        @variables x(t) = 1
        @mtkbuild sys = ODESystem([D(x) ~ x, sspace1.u ~ sspace2.x, sspace2.u ~ sspace1.y],
            t; systems = [sspace1, sspace2])
        prob = ODEProblem(sys, [sspace1.x => 1.0, sspace2.x => 1.0], (0.0, 1.0))
        return sys, prob
    end
    @named sspace1 = StateSpace()
    @named sspace2 = StateSpace()
    truesys, trueprob = build_looped_sspace(sspace1, sspace2)
    truesol = solve(trueprob, Rodas5P(), reltol = 1e-8)
    @test SciMLBase.successful_retcode(truesol)

    @testset "v3, ME" begin
        fmu = loadFMU(joinpath(FMU_DIR, "StateSpace.fmu"); type = :ME)
        @named sspace1 = MTK.FMIComponent(Val(3); fmu, type = :ME)
        @named sspace2 = MTK.FMIComponent(Val(3); fmu, type = :ME)
        sys, prob = build_looped_sspace(sspace1, sspace2)
        sol = solve(prob, Rodas5P(autodiff = false); reltol = 1e-8)
        @test SciMLBase.successful_retcode(sol)
        @test truesol(sol.t;
            idxs = [truesys.sspace1.x, truesys.sspace2.x]).u≈sol(
            sol.t; idxs = [sys.sspace1.x, sys.sspace2.x]).u rtol=1e-7
    end

    @testset "v3, CS" begin
        fmu = loadFMU(joinpath(FMU_DIR, "StateSpace.fmu"); type = :CS)
        @named sspace1 = MTK.FMIComponent(
            Val(3); fmu, type = :CS, communication_step_size = 1e-4)
        @named sspace2 = MTK.FMIComponent(
            Val(3); fmu, type = :CS, communication_step_size = 1e-4)
        sys, prob = build_looped_sspace(sspace1, sspace2)
        sol = solve(prob, Rodas5P(autodiff = false); reltol = 1e-8)
        @test SciMLBase.successful_retcode(sol)
        @test truesol(sol.t;
            idxs = [truesys.sspace1.x, truesys.sspace2.x]).u≈sol(
            sol.t; idxs = [sys.sspace1.x, sys.sspace2.x]).u rtol=1e-2
    end
end
