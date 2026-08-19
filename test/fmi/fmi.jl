using ModelingToolkit, FMI, FMIZoo, OrdinaryDiffEq, NonlinearSolve, SciMLBase
using ModelingToolkit: t_nounits as t, D_nounits as D
using DiffEqBase: BrownFullBasicInit
import ModelingToolkit as MTK

const FMU_DIR = joinpath(@__DIR__, "fmus")

@testset "Standalone pendulum model" begin
    fmu = loadFMU("SpringPendulum1D", "Dymola", "2022x"; type = :ME)
    truesol = FMI.simulate(
        fmu, (0.0, 8.0); saveat = 0.0:0.1:8.0, recordValues = ["mass.s", "mass.v"],
        reltol = 1.0e-8, abstol = 1.0e-8
    )

    function test_no_inputs_outputs(sys)
        for var in unknowns(sys)
            @test !MTK.isinput(var)
            @test !MTK.isoutput(var)
        end
    end
    @testset "v2, ME" begin
        fmu = loadFMU("SpringPendulum1D", "Dymola", "2022x"; type = :ME)
        @mtkcompile sys = MTK.FMIComponent(Val(2); fmu, type = :ME)
        test_no_inputs_outputs(sys)
        prob = ODEProblem{true, SciMLBase.FullSpecialize}(
            sys, [sys.mass__s => 0.5, sys.mass__v => 0.0], (0.0, 8.0)
        )
        sol = solve(prob, Tsit5(); reltol = 1.0e-8, abstol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)

        @test sol(
            0.0:0.1:8.0;
            idxs = [sys.mass__s, sys.mass__v]
        ).u ≈ collect.(truesol.values.saveval) atol = 1.0e-4
        # repeated solve works
        @test_nowarn solve(prob, Tsit5())
    end
    @testset "v2, CS" begin
        fmu = loadFMU("SpringPendulum1D", "Dymola", "2022x"; type = :CS)
        @named inner = MTK.FMIComponent(
            Val(2); fmu, communication_step_size = 1.0e-5, type = :CS
        )
        @variables x(t) = 1.0
        @mtkcompile sys = System([D(x) ~ x], t; systems = [inner])
        test_no_inputs_outputs(sys)

        prob = ODEProblem{true, SciMLBase.FullSpecialize}(
            sys, [sys.inner.mass__s => 0.5, sys.inner.mass__v => 0.0], (0.0, 8.0)
        )
        sol = solve(prob, Tsit5(); reltol = 1.0e-8, abstol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)

        @test sol(
            0.0:0.1:8.0;
            idxs = [sys.inner.mass__s, sys.inner.mass__v]
        ).u ≈ collect.(truesol.values.saveval) rtol = 1.0e-2
    end

    fmu = loadFMU("SpringPendulum1D", "Dymola", "2023x", "3.0"; type = :ME)
    truesol = FMI.simulate(
        fmu, (0.0, 8.0); solver = Tsit5(), saveat = 0.0:0.1:8.0, recordValues = ["mass.s", "mass.v"], tstops = collect(0.0:0.1:8.0)
    )
    @testset "v3, ME" begin
        fmu = loadFMU("SpringPendulum1D", "Dymola", "2023x", "3.0"; type = :ME)
        @mtkcompile sys = MTK.FMIComponent(Val(3); fmu, type = :ME)
        test_no_inputs_outputs(sys)
        prob = ODEProblem{true, SciMLBase.FullSpecialize}(
            sys, [sys.mass__s => 0.5, sys.mass__v => 0.0], (0.0, 8.0)
        )
        sol = solve(prob, Tsit5(); reltol = 1.0e-8, abstol = 1.0e-8, tstops = collect(0.0:0.1:8.0))
        @test SciMLBase.successful_retcode(sol)

        @test sol(
            0.0:0.1:8.0;
            idxs = [sys.mass__s, sys.mass__v]
        ).u ≈ collect.(truesol.values.saveval) atol = 1.0e-4
        # repeated solve works
        @test_nowarn solve(prob, Tsit5())
    end
    @testset "v3, CS" begin
        fmu = loadFMU("SpringPendulum1D", "Dymola", "2023x", "3.0"; type = :CS)
        @named inner = MTK.FMIComponent(
            Val(3); fmu, communication_step_size = 1.0e-5, type = :CS
        )
        @variables x(t) = 1.0
        @mtkcompile sys = System([D(x) ~ x], t; systems = [inner])
        test_no_inputs_outputs(sys)

        prob = ODEProblem{true, SciMLBase.FullSpecialize}(
            sys, [sys.inner.mass__s => 0.5, sys.inner.mass__v => 0.0], (0.0, 8.0)
        )
        sol = solve(prob, Tsit5(); reltol = 1.0e-8, abstol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)

        @test sol(
            0.0:0.1:8.0;
            idxs = [sys.inner.mass__s, sys.inner.mass__v]
        ).u ≈ collect.(truesol.values.saveval) rtol = 1.0e-2
    end
end

@component function SimpleAdder(; name)
    vars = @variables begin
        a(t)
        b(t)
        c(t)
        out(t)
        out2(t)
    end
    pars = @parameters begin
        value = 1.0
    end
    eqs = [
        out ~ a + b + value
        D(c) ~ out
        out2 ~ 2c
    ]
    System(eqs, t, vars, pars; name)
end

@component function StateSpace(; name)
    vars = @variables begin
        x(t)
        y(t)
        u(t)
    end
    pars = @parameters begin
        A = 1.0
        B = 1.0
        C = 1.0
        _D = 1.0
    end
    eqs = [
        D(x) ~ A * x + B * u
        y ~ C * x + _D * u
    ]
    System(eqs, t, vars, pars; name)
end

@testset "IO Model" begin
    function build_simple_adder(adder)
        @variables a(t) b(t) c(t) [guess = 1.0]
        @mtkcompile sys = System(
            [
                adder.a ~ a, adder.b ~ b, D(a) ~ t,
                D(b) ~ adder.out + adder.c, c^2 ~ adder.out + adder.value,
            ],
            t;
            systems = [adder]
        )
        # c will be solved for by initialization
        # this tests that initialization also works with FMUs
        prob = ODEProblem(
            sys, [sys.adder.c => 2.0, sys.a => 1.0, sys.b => 1.0, sys.adder.value => 2.0],
            (0.0, 1.0)
        )
        return sys, prob
    end

    @named adder = SimpleAdder()
    truesys, trueprob = build_simple_adder(adder)
    truesol = solve(trueprob, abstol = 1.0e-8, reltol = 1.0e-8)
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
        sol = solve(prob, Rodas5P(autodiff = AutoFiniteDiff()), abstol = 1.0e-8, reltol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)

        @test truesol(
            sol.t;
            idxs = [truesys.a, truesys.b, truesys.c, truesys.adder.c]
        ).u ≈ sol[
            [
                sys.a, sys.b, sys.c, sys.adder.c,
            ],
        ] rtol = 1.0e-7
    end
    @testset "v2, CS" begin
        fmu = loadFMU(joinpath(FMU_DIR, "SimpleAdder.fmu"); type = :CS)
        @named adder = MTK.FMIComponent(
            Val(2); fmu, type = :CS, communication_step_size = 1.0e-6,
            reinitializealg = BrownFullBasicInit(abstol = 1.0e-7)
        )
        @test MTK.isinput(adder.a)
        @test MTK.isinput(adder.b)
        @test MTK.isoutput(adder.out)
        @test MTK.isoutput(adder.out2)
        @test !MTK.isinput(adder.c) && !MTK.isoutput(adder.c)

        sys, prob = build_simple_adder(adder)
        sol = solve(prob, Rodas5P(autodiff = AutoFiniteDiff()), abstol = 1.0e-8, reltol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)

        @test truesol(
            sol.t;
            idxs = [truesys.a, truesys.b, truesys.c]
        ).u ≈ sol[[sys.a, sys.b, sys.c]] rtol = 4.0e-2
        # sys.adder.c is a discrete variable
        @test truesol(sol.t; idxs = truesys.adder.c).u ≈ sol(sol.t; idxs = sys.adder.c).u rtol = 4.0e-2
    end

    function build_sspace_model(sspace)
        @variables u(t) = 1.0 x(t) = 1.0 y(t) [guess = 1.0]
        @mtkcompile sys = System(
            [sspace.u ~ u, D(u) ~ t, D(x) ~ sspace.x + sspace.y, y^2 ~ sspace.y + sspace.x], t;
            systems = [sspace]
        )

        prob = ODEProblem(
            sys, [sys.sspace.x => 1.0, sys.sspace.A => 2.0], (0.0, 1.0); use_scc = false
        )
        return sys, prob
    end

    @named sspace = StateSpace()
    truesys, trueprob = build_sspace_model(sspace)
    truesol = solve(trueprob, abstol = 1.0e-8, reltol = 1.0e-8)
    @test SciMLBase.successful_retcode(truesol)

    @testset "v3, ME" begin
        fmu = loadFMU(joinpath(FMU_DIR, "StateSpace.fmu"); type = :ME)
        @named sspace = MTK.FMIComponent(Val(3); fmu, type = :ME)
        @test MTK.isinput(sspace.u)
        @test MTK.isoutput(sspace.y)
        @test !MTK.isinput(sspace.x) && !MTK.isoutput(sspace.x)

        sys, prob = build_sspace_model(sspace)
        sol = solve(prob, Rodas5P(autodiff = AutoFiniteDiff()); abstol = 1.0e-8, reltol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)

        @test truesol(
            sol.t;
            idxs = [truesys.u, truesys.x, truesys.y, truesys.sspace.x]
        ).u ≈ sol[
            [
                sys.u, sys.x, sys.y, sys.sspace.x,
            ],
        ] rtol = 1.0e-7
    end

    @testset "v3, CS" begin
        fmu = loadFMU(joinpath(FMU_DIR, "StateSpace.fmu"); type = :CS)
        @named sspace = MTK.FMIComponent(
            Val(3); fmu, communication_step_size = 1.0e-6, type = :CS,
            reinitializealg = BrownFullBasicInit(abstol = 1.0e-7)
        )
        @test MTK.isinput(sspace.u)
        @test MTK.isoutput(sspace.y)
        @test !MTK.isinput(sspace.x) && !MTK.isoutput(sspace.x)

        sys, prob = build_sspace_model(sspace)
        sol = solve(prob, Rodas5P(autodiff = AutoFiniteDiff()); abstol = 1.0e-8, reltol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)

        @test truesol(
            sol.t; idxs = [truesys.u, truesys.x, truesys.y]
        ).u ≈ sol[[sys.u, sys.x, sys.y]] rtol = 1.0e-2
        @test truesol(sol.t; idxs = truesys.sspace.x).u ≈ sol(sol.t; idxs = sys.sspace.x).u rtol = 1.0e-2
    end
end

@testset "FMUs in a loop" begin
    function build_looped_adders(adder1, adder2)
        @variables x(t) = 1
        @mtkcompile sys = System(
            [
                D(x) ~ x, adder1.a ~ adder2.out2,
                adder2.a ~ adder1.out2, adder1.b ~ 1.0, adder2.b ~ 2.0,
            ],
            t;
            systems = [adder1, adder2]
        )
        prob = ODEProblem(
            sys, [adder1.c => 1.0, adder2.c => 1.0, adder1.a => 2.0],
            (0.0, 1.0); guesses = [adder2.a => 0.0]
        )
        return sys, prob
    end
    @named adder1 = SimpleAdder()
    @named adder2 = SimpleAdder()
    truesys, trueprob = build_looped_adders(adder1, adder2)
    truesol = solve(trueprob, Tsit5(), reltol = 1.0e-8)
    @test SciMLBase.successful_retcode(truesol)

    @testset "v2, ME" begin
        fmu = loadFMU(joinpath(FMU_DIR, "SimpleAdder.fmu"); type = :ME)
        @named adder1 = MTK.FMIComponent(Val(2); fmu, type = :ME)
        @named adder2 = MTK.FMIComponent(Val(2); fmu, type = :ME)
        sys, prob = build_looped_adders(adder1, adder2)
        sol = solve(prob, Rodas5P(autodiff = AutoFiniteDiff()); reltol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)
        @test truesol(
            sol.t;
            idxs = [truesys.adder1.c, truesys.adder2.c]
        ).u ≈ sol(
            sol.t; idxs = [sys.adder1.c, sys.adder2.c]
        ).u rtol = 1.0e-7
    end
    @testset "v2, CS" begin
        fmu = loadFMU(joinpath(FMU_DIR, "SimpleAdder.fmu"); type = :CS)
        @named adder1 = MTK.FMIComponent(
            Val(2); fmu, type = :CS, communication_step_size = 1.0e-5
        )
        @named adder2 = MTK.FMIComponent(
            Val(2); fmu, type = :CS, communication_step_size = 1.0e-5
        )
        sys, prob = build_looped_adders(adder1, adder2)
        sol = solve(
            prob,
            Tsit5();
            reltol = 1.0e-8,
            initializealg = SciMLBase.OverrideInit(nlsolve = FastShortcutNLLSPolyalg(autodiff = AutoFiniteDiff()))
        )
        @test truesol(
            sol.t;
            idxs = [truesys.adder1.c, truesys.adder2.c]
        ).u ≈ sol(
            sol.t; idxs = [sys.adder1.c, sys.adder2.c]
        ).u rtol = 1.0e-3
    end

    @testset "multiDimArray Support" begin
        path_to_FMU = joinpath(FMU_DIR, "SimpleArrayModel.fmu")
        fmu = loadFMU(path_to_FMU)
        @named model = MTK.FMIComponent(Val(2); fmu, type = :ME)
        @test model !== nothing
    end
    function build_looped_sspace(sspace1, sspace2)
        @variables x(t) = 1
        @mtkcompile sys = System(
            [D(x) ~ x, sspace1.u ~ sspace2.x, sspace2.u ~ sspace1.y],
            t; systems = [sspace1, sspace2]
        )
        prob = ODEProblem(sys, [sspace1.x => 1.0, sspace2.x => 1.0], (0.0, 1.0))
        return sys, prob
    end
    @named sspace1 = StateSpace()
    @named sspace2 = StateSpace()
    truesys, trueprob = build_looped_sspace(sspace1, sspace2)
    truesol = solve(trueprob, Rodas5P(), reltol = 1.0e-8)
    @test SciMLBase.successful_retcode(truesol)

    @testset "v3, ME" begin
        fmu = loadFMU(joinpath(FMU_DIR, "StateSpace.fmu"); type = :ME)
        @named sspace1 = MTK.FMIComponent(Val(3); fmu, type = :ME)
        @named sspace2 = MTK.FMIComponent(Val(3); fmu, type = :ME)
        sys, prob = build_looped_sspace(sspace1, sspace2)
        sol = solve(prob, Rodas5P(autodiff = AutoFiniteDiff()); reltol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)
        @test truesol(
            sol.t;
            idxs = [truesys.sspace1.x, truesys.sspace2.x]
        ).u ≈ sol(
            sol.t; idxs = [sys.sspace1.x, sys.sspace2.x]
        ).u rtol = 1.0e-7
    end

    @testset "v3, CS" begin
        fmu = loadFMU(joinpath(FMU_DIR, "StateSpace.fmu"); type = :CS)
        @named sspace1 = MTK.FMIComponent(
            Val(3); fmu, type = :CS, communication_step_size = 1.0e-5
        )
        @named sspace2 = MTK.FMIComponent(
            Val(3); fmu, type = :CS, communication_step_size = 1.0e-5
        )
        sys, prob = build_looped_sspace(sspace1, sspace2)
        sol = solve(prob, Rodas5P(autodiff = AutoFiniteDiff()); reltol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)
        @test truesol(
            sol.t;
            idxs = [truesys.sspace1.x, truesys.sspace2.x]
        ).u ≈ sol(
            sol.t; idxs = [sys.sspace1.x, sys.sspace2.x]
        ).u rtol = 1.0e-2
    end
end

@testset "FMU event metadata" begin
    ext = Base.get_extension(ModelingToolkit, :MTKFMIExt)
    # every parameter of `sys` carrying `FMUEventMetadata`, with its metadata
    function event_metadata(sys)
        return [
            (par, MTK.getmetadata(par, ext.FMUEventMetadata))
                for par in parameters(sys)
                if MTK.hasmetadata(par, ext.FMUEventMetadata)
        ]
    end

    @testset "v2, ME" begin
        fmu = loadFMU("BouncingBall1D", "Dymola", "2023x"; type = :ME)
        @named ball = MTK.FMIComponent(Val(2); fmu, type = :ME)
        @variables x(t) = 1.0
        @mtkcompile sys = System([D(x) ~ -x], t; systems = [ball])
        par, meta = only(event_metadata(sys))
        @test MTK.getname(par) == :ball₊wrapper
        @test meta.fmi_version == 2
        @test meta.interface == :ME
        @test meta.n_event_indicators == 2
        @test meta.can_have_time_events
        @test !meta.cs_has_event_mode
        @test Set(zip(meta.state_names, meta.state_vrs)) ==
            Set([(:mass_s, 0x02000000), (:mass_v, 0x02000001)])
        @test isempty(meta.input_names)
        @test ext.resolve_relative(par, meta.state_names[1]) in MTK.getname.(unknowns(sys))
    end

    @testset "v2, CS" begin
        fmu = loadFMU("BouncingBall1D", "Dymola", "2023x"; type = :CS)
        @named ball = MTK.FMIComponent(
            Val(2); fmu, type = :CS, communication_step_size = 1.0e-3
        )
        @variables x(t) = 1.0
        @mtkcompile sys = System([D(x) ~ -x], t; systems = [ball])
        par, meta = only(event_metadata(sys))
        @test MTK.getname(par) == :ball₊wrapper
        @test meta.fmi_version == 2
        @test meta.interface == :CS
        @test meta.n_event_indicators == 2
        @test !meta.can_have_time_events
        # FMI2 Co-Simulation has no event mode
        @test !meta.cs_has_event_mode
        @test Set(meta.state_names) == Set([:mass_s, :mass_v])
        @test ext.resolve_relative(par, :mass_s) === :ball₊mass_s
    end

    @testset "v3, ME" begin
        fmu = loadFMU("SpringFrictionPendulum1D", "Dymola", "2023x", "3.0"; type = :ME)
        @named pendulum = MTK.FMIComponent(Val(3); fmu, type = :ME)
        @variables x(t) = 1.0
        @mtkcompile sys = System([D(x) ~ -x], t; systems = [pendulum])
        par, meta = only(event_metadata(sys))
        @test MTK.getname(par) == :pendulum₊wrapper
        @test meta.fmi_version == 3
        @test meta.interface == :ME
        @test meta.n_event_indicators == 32
        @test meta.can_have_time_events
        @test !meta.cs_has_event_mode
        @test Set(zip(meta.state_names, meta.state_vrs)) ==
            Set([(:mass__s, 0x02000000), (:mass__v, 0x02000001)])
        @test ext.resolve_relative(par, meta.state_names[1]) in MTK.getname.(unknowns(sys))
    end

    @testset "v3, CS" begin
        fmu = loadFMU("SpringFrictionPendulum1D", "Dymola", "2023x", "3.0"; type = :CS)
        @named pendulum = MTK.FMIComponent(
            Val(3); fmu, type = :CS, communication_step_size = 1.0e-3
        )
        @variables x(t) = 1.0
        @mtkcompile sys = System([D(x) ~ -x], t; systems = [pendulum])
        par, meta = only(event_metadata(sys))
        @test MTK.getname(par) == :pendulum₊wrapper
        @test meta.fmi_version == 3
        @test meta.interface == :CS
        @test meta.n_event_indicators == 32
        @test !meta.can_have_time_events
        @test meta.cs_has_event_mode
        @test Set(meta.state_names) == Set([:mass__s, :mass__v])
        @test ext.resolve_relative(par, :mass__v) === :pendulum₊mass__v
    end

    @testset "two FMUs in one system" begin
        fmu = loadFMU("BouncingBall1D", "Dymola", "2023x"; type = :ME)
        @named ball1 = MTK.FMIComponent(Val(2); fmu, type = :ME)
        @named ball2 = MTK.FMIComponent(Val(2); fmu, type = :ME)
        @variables x(t) = 1.0
        @mtkcompile sys = System([D(x) ~ -x], t; systems = [ball1, ball2])
        carriers = event_metadata(sys)
        @test length(carriers) == 2
        @test Set(MTK.getname(par) for (par, _) in carriers) ==
            Set([:ball1₊wrapper, :ball2₊wrapper])
        # one metadata object per component, not a shared one
        @test carriers[1][2] !== carriers[2][2]
        unknown_names = MTK.getname.(unknowns(sys))
        for (par, meta) in carriers
            @test meta.n_event_indicators == 2
            @test Set(meta.state_names) == Set([:mass_s, :mass_v])
            for name in meta.state_names
                @test ext.resolve_relative(par, name) in unknown_names
            end
        end
    end

    @testset "inputs, no event indicators, un-namespaced wrapper" begin
        fmu = loadFMU(joinpath(FMU_DIR, "StateSpace.fmu"); type = :ME)
        @named sspace = MTK.FMIComponent(Val(3); fmu, type = :ME)
        par, meta = only(event_metadata(sspace))
        @test MTK.getname(par) == :wrapper
        @test meta.n_event_indicators == 0
        @test meta.state_names == [:x]
        @test meta.state_vrs == [0x0000000b]
        @test meta.input_names == [:u]
        @test ext.resolve_relative(par, :x) === :x
    end
end
