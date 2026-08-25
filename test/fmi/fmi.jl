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
    # the variables of `component` that are continuous states of its FMU
    function state_variables(component, meta)
        vars = [var for var in unknowns(component) if MTK.getname(var) in meta.state_names]
        @test length(vars) == length(meta.state_names)
        return vars
    end
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
        @test Set(zip(meta.state_names, meta.state_value_references)) ==
            Set([(:mass_s, 0x02000000), (:mass_v, 0x02000001)])
        @test isempty(meta.input_names)
        @test ext.resolve_relative(par, meta.state_names[1]) in MTK.getname.(unknowns(sys))
        # a Model Exchange FMU that declares event indicators can have events, and an event
        # writes the values it produces for the continuous states into the integrator's `u`,
        # so the import keeps them out of simplification
        @test all(MTK.isirreducible, state_variables(ball, meta))
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
        # a CoSimulation FMU integrates its own states, which never enter the integrator's
        # `u`, so nothing keeps them out of simplification
        @test !any(MTK.isirreducible, state_variables(ball, meta))
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
        @test Set(zip(meta.state_names, meta.state_value_references)) ==
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
        @test meta.state_value_references == [0x0000000b]
        @test meta.input_names == [:u]
        @test ext.resolve_relative(par, :x) === :x
        # an FMU with no event indicators has no state event to write back, so its states are
        # left to simplification and it falls back to being notified of completed steps
        @test !any(MTK.isirreducible, state_variables(sspace, meta))
    end
end

@testset "state-valued outputs" begin
    # Reference-FMUs `BouncingBall` declares `h` and `v` both as continuous states and as
    # outputs. The analytic solution below holds until the first bounce at ≈0.45s.
    analytic_h(tval) = 1.0 - 9.81 / 2 * tval^2
    analytic_v(tval) = -9.81 * tval

    function build_ball_model(ball)
        @variables x(t) = 1.0
        @mtkcompile sys = System([D(x) ~ -x], t; systems = [ball])
        prob = ODEProblem(sys, [ball.h => 1.0, ball.v => 0.0], (0.0, 0.2))
        return sys, prob
    end

    @testset "v$Ver, ME" for Ver in (2, 3)
        fmu = loadFMU(joinpath(FMU_DIR, "BouncingBall$Ver.fmu"); type = :ME)
        @named ball = MTK.FMIComponent(Val(Ver); fmu, type = :ME)
        # the state variable is reused as the output, so it keeps its own causality
        @test !MTK.isoutput(ball.h)
        sys, prob = build_ball_model(ball)
        unknown_names = MTK.getname.(unknowns(sys))
        @test count(isequal(:ball₊h), unknown_names) == 1
        @test count(isequal(:ball₊v), unknown_names) == 1

        sol = solve(prob, Tsit5(); abstol = 1.0e-8, reltol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)
        for tval in (0.05, 0.1, 0.2)
            @test sol(tval; idxs = ball.h) ≈ analytic_h(tval) atol = 1.0e-3
        end
        @test sol(0.2; idxs = ball.v) ≈ analytic_v(0.2) atol = 1.0e-3
    end

    @testset "v$Ver, CS" for Ver in (2, 3)
        fmu = loadFMU(joinpath(FMU_DIR, "BouncingBall$Ver.fmu"); type = :CS)
        @named ball = MTK.FMIComponent(
            Val(Ver); fmu, type = :CS, communication_step_size = 1.0e-3
        )
        sys, prob = build_ball_model(ball)
        sol = solve(prob, Tsit5(); abstol = 1.0e-8, reltol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)
        # the FMU integrates itself with a fixed step, hence the looser tolerance
        @test sol(0.2; idxs = ball.h) ≈ analytic_h(0.2) atol = 1.0e-2
        @test sol(0.2; idxs = ball.v) ≈ analytic_v(0.2) atol = 1.0e-2
    end
end

@testset "FMU event iteration" begin
    ext = Base.get_extension(ModelingToolkit, :MTKFMIExt)

    # the instance wrapper is the default value of the callable parameter `FMIComponent`
    # creates for it
    function ball_wrapper(Ver; kwargs...)
        fmu = loadFMU(joinpath(FMU_DIR, "BouncingBall$Ver.fmu"); type = :ME)
        @named ball = MTK.FMIComponent(Val(Ver); fmu, type = :ME, kwargs...)
        return MTK.getdefault(
            only(p for p in parameters(ball) if MTK.getname(p) == :wrapper)
        )
    end

    # an FMU enters Event Mode when it leaves initialization mode, which is where
    # `do_fmu_event_iteration!` expects to be called. Empty `params` leaves the FMU's own
    # start values in place.
    function enter_initial_event_mode!(wrapper::ext.FMI2InstanceWrapper)
        ext.get_instance_common!(wrapper, Float64[], Float64[], 0.0)
        ext.FMI.fmi2ExitInitializationMode(wrapper.instance)
        return wrapper
    end
    function enter_initial_event_mode!(wrapper::ext.FMI3InstanceWrapper)
        wrapper.instance = ext.FMI.fmi3InstantiateModelExchange!(wrapper.fmu)
        ext.get_instance_common!(wrapper, Float64[], Float64[], 0.0)
        ext.FMI.fmi3ExitInitializationMode(wrapper.instance)
        return wrapper
    end

    @testset "v$Ver, ME" for Ver in (2, 3)
        wrapper = ball_wrapper(Ver)
        @test wrapper.next_event_time === nothing
        @test length(wrapper.event_states_buffer) == length(wrapper.state_value_references)

        enter_initial_event_mode!(wrapper)
        result = ext.do_fmu_event_iteration!(wrapper)
        @test !result.terminate
        # `BouncingBall` starts at h = 1 with no event pending
        @test !result.values_changed
        @test !result.nominals_changed
        @test result.next_event_time === nothing
        ext.reset_instance!(wrapper)
        @test wrapper.instance === nothing

        # the same iteration as part of instance initialization
        ext.get_instance_ME!(wrapper, Float64[], Float64[], 0.0)
        @test wrapper.instance !== nothing
        @test wrapper.next_event_time === nothing
        @test ext.partiallyCompleteIntegratorStep(wrapper) ==
            (; enter_event_mode = false, terminate = false)
        ext.reset_instance!(wrapper)
    end

    @testset "the iteration cap is per FMU, v$Ver" for Ver in (2, 3)
        # the cap lives on the instance wrapper, so every FMU carries its own. No FMU in the
        # test suite fails to converge, so the kwarg is checked where it lands.
        @test ball_wrapper(Ver).max_event_iterations == 100
        @test ball_wrapper(Ver; max_event_iterations = 7).max_event_iterations == 7
        # a cap below one runs no iteration at all, which would report the FMU as
        # non-convergent instead of the cap as nonsense
        @test_throws ArgumentError ball_wrapper(Ver; max_event_iterations = 0)
    end

    @testset "status handling" begin
        @test !ext.fmu_status_is_error(Val(2), nothing, :fmi2SetReal)
        @test !ext.fmu_status_is_error(Val(2), FMI.fmi2StatusOK, :fmi2SetReal)
        @test !ext.fmu_status_is_error(Val(3), FMI.fmi3StatusOK, :fmi3SetFloat64)
        # a warning is logged but is not a failure
        @test @test_logs (:warn,) !ext.fmu_status_is_error(
            Val(2), FMI.fmi2StatusWarning, :fmi2SetReal
        )
        @test @test_logs (:warn,) !ext.fmu_status_is_error(
            Val(3), FMI.fmi3StatusWarning, :fmi3SetFloat64
        )
        @test ext.fmu_status_is_error(Val(2), FMI.fmi2StatusError, :fmi2SetReal)
        @test ext.fmu_status_is_error(Val(2), FMI.fmi2StatusFatal, :fmi2SetReal)
        @test ext.fmu_status_is_error(Val(3), FMI.fmi3StatusError, :fmi3SetFloat64)
        @test ext.fmu_status_is_error(Val(3), FMI.fmi3StatusFatal, :fmi3SetFloat64)

        # the cleanup branch has to use the functions matching the FMI version of the
        # checked call
        leaves(e) = e isa Expr ? mapreduce(leaves, vcat, e.args; init = Any[]) : Any[e]
        v2_expansion = leaves(
            macroexpand(
                ext, :(@statuscheck FMI.fmi2ExitInitializationMode(wrapper.instance))
            )
        )
        @test ext.FMI.fmi2Terminate in v2_expansion
        @test !(ext.FMI.fmi3Terminate in v2_expansion)
        v3_expansion = leaves(
            macroexpand(
                ext, :(@statuscheck FMI.fmi3ExitInitializationMode(wrapper.instance))
            )
        )
        @test ext.FMI.fmi3Terminate in v3_expansion
        @test !(ext.FMI.fmi2Terminate in v3_expansion)
    end
end

@testset "FMU state events" begin
    ext = Base.get_extension(ModelingToolkit, :MTKFMIExt)

    # `≈` on vectors compares 2-norms, which conflates a per-point tolerance with the number
    # of points
    max_error(a, b) = maximum(abs, a .- b)

    # a Reference-FMU ships its own reference output at this path inside the archive
    function reference_output(fmu, model_name)
        path = joinpath(
            fmu.path, "extra", "org.fmi-standard.fmi-ls-ref", "$(model_name)_out.csv"
        )
        lines = collect(eachline(path))
        columns = Symbol.(split(lines[1], ','))
        rows = [parse.(Float64, split(line, ',')) for line in lines[2:end]]
        return columns, permutedims(reduce(hcat, rows))
    end

    @testset "BouncingBall v$Ver, ME" for Ver in (2, 3)
        fmu = loadFMU(joinpath(FMU_DIR, "BouncingBall$Ver.fmu"); type = :ME)
        @named ball = MTK.FMIComponent(Val(Ver); fmu, type = :ME)
        @variables x(t) = 1.0
        @mtkcompile sys = System([D(x) ~ -x], t; systems = [ball])
        prob = ODEProblem(sys, [ball.h => 1.0, ball.v => 0.0], (0.0, 3.0))
        sol = solve(prob, Tsit5(); abstol = 1.0e-8, reltol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)

        heights = sol[ball.h]
        velocities = sol[ball.v]
        # `RightRootFind` localizes the event past the crossing, so `h` may dip below zero by
        # the root-finding tolerance, but never visibly
        @test minimum(heights) >= -1.0e-8
        # `save_positions = (true, true)` duplicates the time of every event
        bounces = findall(i -> sol.t[i] == sol.t[i + 1], 1:(length(sol.t) - 1))
        @test length(bounces) >= 2
        for i in bounces
            # the FMU stops the ball once a bounce would leave it slower than `v_min`, so the
            # last bounce brings `v` to zero rather than to a positive value
            @test velocities[i] < 0 <= velocities[i + 1]
        end

        columns, reference = reference_output(fmu, "BouncingBall")
        @test columns == [:time, :h, :v]
        heights_at_reference = sol(reference[:, 1]; idxs = ball.h).u
        # the reference output is a fixed-step explicit Euler run (`FIXED_SOLVER_STEP = 1e-3`)
        # whose bounce times drift, so the whole-trace bound absorbs the reference's own error
        @test max_error(heights_at_reference, reference[:, 2]) <= 3.0e-2
        # up to the first bounce the reference is a plain free fall and is accurate
        free_fall = reference[:, 1] .<= 0.4
        @test max_error(heights_at_reference[free_fall], reference[free_fall, 2]) <= 1.0e-2
    end

    @testset "BouncingBall v3, CS: event mode" begin
        fmu = loadFMU(joinpath(FMU_DIR, "BouncingBall3.fmu"); type = :CS)
        @test ext.cs_event_mode_supported(fmu)
        # the FMU's own fixed solver step, so every one of its internal steps ends at a
        # communication point and every state event is reported to the importer
        @named ball = MTK.FMIComponent(
            Val(3); fmu, type = :CS, communication_step_size = 1.0e-3
        )
        wrapper = MTK.getdefault(
            only(p for p in parameters(ball) if MTK.getname(p) == :wrapper)
        )
        # an instance using event mode leaves initialization mode in Event Mode; the initial
        # event iteration and `fmi3EnterStepMode` are what put it in Step Mode
        ext.get_instance_CS!(wrapper, Float64[], Float64[], Float64[], 0.0, 3.0)
        @test wrapper.instance.state == ext.FMI.fmi3InstanceStateStepMode
        @test wrapper.next_event_time === nothing
        ext.reset_instance!(wrapper)

        @variables x(t) = 1.0
        @mtkcompile sys = System([D(x) ~ -x], t; systems = [ball])
        par = only(p for p in parameters(sys) if MTK.hasmetadata(p, ext.FMUEventMetadata))
        @test MTK.getmetadata(par, ext.FMUEventMetadata).cs_has_event_mode
        # event mode lives in the periodic communication callback, so a CoSimulation FMU
        # needs no callback construction hook
        @test !MTK.hasmetadata(par, CallbackConstructionHook)

        prob = ODEProblem(sys, [ball.h => 1.0, ball.v => 0.0], (0.0, 3.0))
        sol = solve(prob, Tsit5(); abstol = 1.0e-8, reltol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)

        columns, reference = reference_output(fmu, "BouncingBall")
        @test columns == [:time, :h, :v]
        heights = sol(reference[:, 1]; idxs = ball.h).u
        velocities = sol(reference[:, 1]; idxs = ball.v).u
        # the reference output is this FMU's own fixed-step trace at the communication points,
        # which matching solver and communication grids reproduce to round-off
        @test max_error(heights, reference[:, 2]) <= 1.0e-6
        @test max_error(velocities, reference[:, 3]) <= 1.0e-6
        @test minimum(heights) >= 0
        # every event the FMU reports has to be handled for the ball to bounce at all: an
        # event left unhandled is lost, and the ball falls through the floor
        bounces(vs) = count(i -> vs[i] < 0 <= vs[i + 1], 1:(length(vs) - 1))
        @test bounces(reference[:, 3]) >= 4
        @test bounces(velocities) == bounces(reference[:, 3])
    end

    @testset "Dahlquist v3, CS: no event mode" begin
        fmu = loadFMU(joinpath(FMU_DIR, "Dahlquist3.fmu"); type = :CS)
        # this FMU declares no event mode, so it is held out of the event-mode path and its
        # trace is unaffected: an event iteration run on it would be rejected as illegal
        @test !ext.cs_event_mode_supported(fmu)
        columns, reference = reference_output(fmu, "Dahlquist")
        @test columns == [:time, :x]
        @named dahlquist = MTK.FMIComponent(
            Val(3); fmu, type = :CS, communication_step_size = 1.0e-1
        )
        @variables y(t) = 1.0
        @mtkcompile sys = System([D(y) ~ -y], t; systems = [dahlquist])
        prob = ODEProblem(sys, [dahlquist.x => 1.0], (0.0, reference[end, 1]))
        sol = solve(prob, Tsit5(); abstol = 1.0e-8, reltol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)
        # the communication step matches the FMU's own solver step and the reference output
        # interval, so the trace is reproduced to round-off. The final communication step of
        # the solve is not saved, hence the last sample is excluded.
        got = sol(reference[1:(end - 1), 1]; idxs = dahlquist.x).u
        @test max_error(got, reference[1:(end - 1), 2]) <= 1.0e-12
    end

    @testset "VanDerPol v2, ME" begin
        fmu = loadFMU(joinpath(FMU_DIR, "VanDerPol2.fmu"); type = :ME)
        @named vdp = MTK.FMIComponent(Val(2); fmu, type = :ME)
        @variables x(t) = 1.0
        @mtkcompile sys = System([D(x) ~ -x], t; systems = [vdp])
        par, meta = only(
            (p, MTK.getmetadata(p, ext.FMUEventMetadata)) for p in parameters(sys)
                if MTK.hasmetadata(p, ext.FMUEventMetadata)
        )
        @test meta.n_event_indicators == 0
        # an FMU with no event indicators has no state-event callback, but every Model
        # Exchange FMU gets the per-accepted-step callback
        @test MTK.hasmetadata(par, CallbackConstructionHook)
        @test only(MTK.getmetadata(par, CallbackConstructionHook)(sys, par)) isa
            SciMLBase.DiscreteCallback

        columns, reference = reference_output(fmu, "VanDerPol")
        @test columns == [:time, :x0, :x1]
        prob = ODEProblem(sys, [vdp.x0 => 2.0, vdp.x1 => 0.0], (0.0, reference[end, 1]))
        sol = solve(prob, Tsit5(); abstol = 1.0e-10, reltol = 1.0e-10)
        @test SciMLBase.successful_retcode(sol)
        x0_at_reference = sol(reference[:, 1]; idxs = vdp.x0).u
        x1_at_reference = sol(reference[:, 1]; idxs = vdp.x1).u
        # the reference output is a 1e-2 Euler step over 20 s whose phase error accumulates
        # across the oscillator's ~6.7 s periods, so the whole-trace bounds absorb the
        # reference's own drift
        @test max_error(x0_at_reference, reference[:, 2]) <= 0.45
        @test max_error(x1_at_reference, reference[:, 3]) <= 0.8
        # over the first 1 s the reference has not drifted yet
        early = reference[:, 1] .<= 1.0
        @test max_error(x0_at_reference[early], reference[early, 2]) <= 1.0e-2
        @test max_error(x1_at_reference[early], reference[early, 3]) <= 1.0e-2
    end

    @testset "Dahlquist v3, ME" begin
        fmu = loadFMU(joinpath(FMU_DIR, "Dahlquist3.fmu"); type = :ME)
        @named dahlquist = MTK.FMIComponent(Val(3); fmu, type = :ME)
        @variables x(t) = 1.0
        @mtkcompile sys = System([D(x) ~ -x], t; systems = [dahlquist])
        columns, reference = reference_output(fmu, "Dahlquist")
        @test columns == [:time, :x]
        prob = ODEProblem(sys, [dahlquist.x => 1.0], (0.0, reference[end, 1]))
        sol = solve(prob, Tsit5(); abstol = 1.0e-10, reltol = 1.0e-10)
        @test SciMLBase.successful_retcode(sol)
        got = sol(reference[:, 1]; idxs = dahlquist.x).u
        # `Dahlquist` is `x' = -k*x` with `k = 1`
        @test max_error(got, exp.(-reference[:, 1])) <= 1.0e-6
        # the reference output uses a 0.1 Euler step
        @test max_error(got, reference[:, 2]) <= 2.5e-2
    end

    @testset "Stair v2, ME: unsupported time events" begin
        fmu = loadFMU(joinpath(FMU_DIR, "Stair2.fmu"); type = :ME)
        @variables x(t) = 1.0
        # `Stair` has no continuous states, so consuming its output is what forces the FMU to
        # be instantiated during the solve
        @named stair = MTK.FMIComponent(Val(2); fmu, type = :ME)
        @mtkcompile sys = System([D(x) ~ -x + stair.counter], t; systems = [stair])
        @test_throws "ignore_time_events" solve(
            ODEProblem(sys, [], (0.0, 3.0)), Tsit5()
        )

        # `Stair`'s only output is an FMI Integer, which the Model Exchange functor cannot
        # read, so the ignoring path is exercised on the instance wrapper directly
        @named lenient_stair = MTK.FMIComponent(
            Val(2); fmu, type = :ME, ignore_time_events = true
        )
        wrapper = MTK.getdefault(
            only(p for p in parameters(lenient_stair) if MTK.getname(p) == :wrapper)
        )
        logger = Test.TestLogger(; min_level = Base.CoreLogging.Warn)
        Base.CoreLogging.with_logger(logger) do
            # only the first call instantiates and runs the initial event iteration, so only
            # it may warn
            ext.get_instance_ME!(wrapper, Float64[], Float64[], 0.0)
            ext.get_instance_ME!(wrapper, Float64[], Float64[], 0.0)
        end
        @test count(
            record -> occursin("time event", string(record.message)), logger.logs
        ) == 1
        @test wrapper.next_event_time == 1.0
        ext.reset_instance!(wrapper)

        # the guard itself, driven directly: instantiating once is what makes the test above
        # see one warning, and the affect-side call is unreachable for an FMU with no event
        # indicators
        wrapper.time_event_warned = false
        logger = Test.TestLogger(; min_level = Base.CoreLogging.Warn)
        Base.CoreLogging.with_logger(logger) do
            ext.handle_fmu_time_event!(wrapper, 1.0)
            ext.handle_fmu_time_event!(wrapper, 2.0)
        end
        @test count(
            record -> occursin("time event", string(record.message)), logger.logs
        ) == 1
        @test wrapper.time_event_warned

        strict_wrapper = MTK.getdefault(
            only(p for p in parameters(stair) if MTK.getname(p) == :wrapper)
        )
        @test_throws "ignore_time_events" ext.handle_fmu_time_event!(strict_wrapper, 1.0)
        @test ext.handle_fmu_time_event!(strict_wrapper, nothing) === nothing
    end

    @testset "Stair v3, ME: unsupported time events" begin
        fmu = loadFMU(joinpath(FMU_DIR, "Stair3.fmu"); type = :ME)
        # `Stair`'s only output is an FMI Integer, which the Model Exchange functor cannot
        # read, so both policies are exercised on the instance wrapper directly
        @named stair = MTK.FMIComponent(Val(3); fmu, type = :ME)
        strict_wrapper = MTK.getdefault(
            only(p for p in parameters(stair) if MTK.getname(p) == :wrapper)
        )
        @test_throws "ignore_time_events" ext.get_instance_ME!(
            strict_wrapper, Float64[], Float64[], 0.0
        )

        @named lenient_stair = MTK.FMIComponent(
            Val(3); fmu, type = :ME, ignore_time_events = true
        )
        wrapper = MTK.getdefault(
            only(p for p in parameters(lenient_stair) if MTK.getname(p) == :wrapper)
        )
        logger = Test.TestLogger(; min_level = Base.CoreLogging.Warn)
        Base.CoreLogging.with_logger(logger) do
            # only the first call instantiates and runs the initial event iteration, so only
            # it may warn
            ext.get_instance_ME!(wrapper, Float64[], Float64[], 0.0)
            ext.get_instance_ME!(wrapper, Float64[], Float64[], 0.0)
        end
        @test count(
            record -> occursin("time event", string(record.message)), logger.logs
        ) == 1
        # `Stair` steps its counter at every integer time, and the v3 event iteration reports
        # the first of those through the `event_time_defined`/`event_time` pair
        @test wrapper.next_event_time == 1.0
        ext.reset_instance!(wrapper)
    end
end

@testset "FMU step events" begin
    ext = Base.get_extension(ModelingToolkit, :MTKFMIExt)

    @testset "the ME functor reads outputs before derivatives" begin
        fmu = loadFMU(joinpath(FMU_DIR, "SimpleAdder.fmu"); type = :ME)
        @named adder = MTK.FMIComponent(Val(2); fmu, type = :ME)
        wrapper = MTK.getdefault(
            only(p for p in parameters(adder) if MTK.getname(p) == :wrapper)
        )
        # `SimpleAdder` is `out ~ a + b + value, D(c) ~ out, out2 ~ 2c`, and the functor
        # returns `[D(c), out2, out]`. `out2` is an output no derivative depends on, which an
        # OpenModelica FMU evaluates only when a variable is read, so reading the derivatives
        # first returns the `out2` of the previous evaluation.
        @test wrapper([1.0], [2.0, 1.0], [1.0], 0.0) ≈ [4.0, 2.0, 4.0]
        @test wrapper([5.0], [2.0, 1.0], [1.0], 0.0) ≈ [4.0, 10.0, 4.0]
        @test wrapper([9.0], [4.0, 1.0], [1.0], 0.25) ≈ [6.0, 18.0, 6.0]
        ext.reset_instance!(wrapper)
    end

    @testset "per-accepted-step notification, BouncingBall v$Ver, ME" for Ver in (2, 3)
        fmu = loadFMU(joinpath(FMU_DIR, "BouncingBall$Ver.fmu"); type = :ME)
        @named ball = MTK.FMIComponent(Val(Ver); fmu, type = :ME)
        @variables x(t) = 1.0
        @mtkcompile sys = System([D(x) ~ -x], t; systems = [ball])
        prob = ODEProblem(sys, [ball.h => 1.0, ball.v => 0.0], (0.0, 3.0))
        sol = solve(prob, Tsit5(); abstol = 1.0e-8, reltol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)
        velocities = sol[ball.v]
        duplicated_times = findall(i -> sol.t[i] == sol.t[i + 1], 1:(length(sol.t) - 1))
        # an accepted step is not a discontinuity, so the notification saves nothing of its
        # own: every duplicated time in the solution is an event, and the FMU bounced the
        # ball at each of them
        @test length(duplicated_times) >= 2
        @test all(i -> velocities[i] < 0 <= velocities[i + 1], duplicated_times)
        @test minimum(sol[ball.h]) >= -1.0e-8
        # the rest of the system is untouched by the notification: `x` decays from 1 for the
        # whole solve
        @test sol[x][end] ≈ exp(-3.0) rtol = 1.0e-5
    end

    @testset "the notification stays out of the FMU where FMI forbids it, v$Ver" for Ver in
        (2, 3)
        # no FMU in the suite requests termination, and an illegal call is refused by the FMU
        # rather than observable in a trajectory, so these two paths are driven directly
        fmu = loadFMU(joinpath(FMU_DIR, "BouncingBall$Ver.fmu"); type = :ME)
        @named ball = MTK.FMIComponent(Val(Ver); fmu, type = :ME)
        @variables x(t) = 1.0
        @mtkcompile sys = System([D(x) ~ -x], t; systems = [ball])
        par = only(p for p in parameters(sys) if MTK.hasmetadata(p, ext.FMUEventMetadata))
        callbacks = MTK.getmetadata(par, CallbackConstructionHook)(sys, par)
        # a reinitialization algorithm leaves no trace in the trajectory of an FMU whose
        # compiled system is a plain ODE, so the resolved default is checked where it lands
        @test all(cb -> cb.initializealg isa SciMLBase.CheckInit, callbacks)

        step_cb = callbacks[2]
        # the FMU is told about the `t0` step through `initialize`, which is what FMI.jl's
        # `func_start` does. No fixture has an event at `t0`, so no trajectory shows it.
        @test step_cb.initialize !== SciMLBase.INITIALIZE_DEFAULT
        prob = ODEProblem(sys, [ball.h => 1.0, ball.v => 0.0], (0.0, 3.0))
        integrator = init(prob, Tsit5(); abstol = 1.0e-8, reltol = 1.0e-8)
        wrapper = ext.SII.getp(sys, par)(integrator)

        # an FMU with no instance has no completed step to report
        ext.reset_instance!(wrapper)
        @test step_cb.affect!(integrator) === nothing
        @test wrapper.instance === nothing

        # a terminating event leaves the FMU in Event Mode, where the calls the notification
        # makes are illegal, and SciML applies discrete callbacks after a continuous one
        # terminated the solve
        ext.get_instance_ME!(wrapper, Float64[], Float64[], integrator.t)
        SciMLBase.terminate!(integrator)
        ext.enter_fmu_event_mode!(wrapper, nothing)
        @test step_cb.affect!(integrator) === nothing
        # the FMU was left alone, so it is still in Event Mode
        @test ext.leave_fmu_event_mode!(wrapper, false) === nothing
        ext.reset_instance!(wrapper)
    end

    @testset "step event iteration, BouncingBall v$Ver" for Ver in (2, 3)
        fmu = loadFMU(joinpath(FMU_DIR, "BouncingBall$Ver.fmu"); type = :ME)
        @named ball = MTK.FMIComponent(Val(Ver); fmu, type = :ME)
        wrapper = MTK.getdefault(
            only(p for p in parameters(ball) if MTK.getname(p) == :wrapper)
        )
        ext.get_instance_ME!(wrapper, Float64[], Float64[], 0.0)
        # the FMU keeps its own start values throughout initialization, so the event
        # iteration performed there cannot see the state MTK integrates from
        @test wrapper.event_states_buffer == [0.0, 0.0]

        ext.force_set_fmu_state!(wrapper, [-1.0, -1.0], Float64[], 0.0)
        # a step event has no triggered event indicator
        ext.enter_fmu_event_mode!(wrapper, nothing)
        result = ext.do_fmu_event_iteration!(wrapper)
        @test result.values_changed
        @test !result.terminate
        states = ext.leave_fmu_event_mode!(wrapper, result.values_changed)
        @test states === wrapper.event_states_buffer
        # `BouncingBall` reflects the ball: `h` just above zero, `v` reversed and damped
        @test states[1] > 0
        @test states[2] ≈ 0.7
        ext.reset_instance!(wrapper)
    end

    @testset "event reinitialization on a DAE coupling, BouncingBall v$Ver" for Ver in
        (2, 3)
        fmu = loadFMU(joinpath(FMU_DIR, "BouncingBall$Ver.fmu"); type = :ME)
        # `v` enters the algebraic equation below, so the reversal the FMU writes back at a
        # bounce leaves that equation unsatisfied: the default `CheckInit` errors at the first
        # bounce, while `BrownFullBasicInit` re-solves `z` and the solve runs through
        @named ball = MTK.FMIComponent(
            Val(Ver); fmu, type = :ME, reinitializealg = BrownFullBasicInit()
        )
        # a cubic has no symbolic solution for `z`, so `z` survives simplification as an
        # algebraic unknown and the compiled system gets a singular mass matrix
        @variables z(t) [guess = 1.0]
        @mtkcompile sys = System([z^3 + z ~ ball.v + 2.0], t; systems = [ball])
        prob = ODEProblem(sys, [ball.h => 1.0, ball.v => 0.0], (0.0, 3.0))
        # guard that the system stays a DAE, so the callback reinitialization path is exercised
        mass_matrix = prob.f.mass_matrix
        @test any(i -> iszero(mass_matrix[i, i]), axes(mass_matrix, 1))
        alg = Rodas5P(autodiff = AutoFiniteDiff())
        integrator = init(prob, alg; abstol = 1.0e-8, reltol = 1.0e-8)
        @test integrator.isdae

        sol = solve(prob, alg; abstol = 1.0e-8, reltol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)
        heights = sol[ball.h]
        velocities = sol[ball.v]
        @test minimum(heights) >= -1.0e-8
        # `save_positions = (true, true)` duplicates the time of every event
        bounces = findall(i -> sol.t[i] == sol.t[i + 1], 1:(length(sol.t) - 1))
        @test length(bounces) >= 2
        # the velocity the FMU wrote back has to survive the reinitialization that follows a
        # discontinuity: `OverrideInit` resets `u` to `u0`, which would leave the ball falling
        # from its initial state after every bounce instead of rebounding
        @test count(i -> velocities[i] < 0 < velocities[i + 1], bounces) >= 2
        # the algebraic constraint holds at every saved point, events included
        @test maximum(abs, sol[z] .^ 3 .+ sol[z] .- velocities .- 2.0) <= 1.0e-6

        # the same coupling on the default algorithm, which only checks the constraint
        @named default_ball = MTK.FMIComponent(Val(Ver); fmu, type = :ME)
        @mtkcompile default_sys = System(
            [z^3 + z ~ default_ball.v + 2.0], t; systems = [default_ball]
        )
        default_prob = ODEProblem(
            default_sys, [default_ball.h => 1.0, default_ball.v => 0.0], (0.0, 3.0)
        )
        err = try
            solve(default_prob, alg; abstol = 1.0e-8, reltol = 1.0e-8)
        catch e
            e
        end
        @test err isa SciMLBase.CheckInitFailureError
        # the residual is the velocity reversal itself, not interpolation error at the root
        @test err.normresid > 1.0
    end

    @testset "a 0-indicator FMU whose state cannot be written back" begin
        fmu = loadFMU(joinpath(FMU_DIR, "VanDerPol2.fmu"); type = :ME)
        @named vdp = MTK.FMIComponent(Val(2); fmu, type = :ME)
        @variables y(t) = 1.0
        # declaring a continuous state of the FMU an input of the compiled system takes it out
        # of the unknowns, so a value an event produced for it would have nowhere to go
        sys = MTK.mtkcompile(
            System([D(y) ~ -y], t; systems = [vdp], name = :outer); inputs = [vdp.x1]
        )
        @test ext.SII.variable_index(sys, vdp.x1) === nothing
        par = only(p for p in parameters(sys) if MTK.hasmetadata(p, ext.FMUEventMetadata))
        # an FMU declaring no event indicators only has to be told about completed steps, so
        # it gets a step callback that pushes nothing rather than a construction-time error
        step_cb = only(MTK.getmetadata(par, CallbackConstructionHook)(sys, par))
        @test step_cb isa SciMLBase.DiscreteCallback
        prob = ODEProblem(sys, [vdp.x0 => 2.0, vdp.x1 => 0.0], (0.0, 1.0))
        sol = solve(prob, Tsit5(); abstol = 1.0e-8, reltol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)
        @test sol[y][end] ≈ exp(-1.0)

        # whether the FMU's inputs can be evaluated is only answerable by building their
        # getter: `is_observed` is true even for a symbol the system does not contain, so a
        # predicate cannot tell the two apart
        @test ext.SII.is_observed(sys, :totally_bogus)
        @test_throws ArgumentError ext.SII.observed(sys, [:totally_bogus])
    end

    @testset "an FMU with event indicators keeps its states writable" begin
        fmu = loadFMU(joinpath(FMU_DIR, "BouncingBall2.fmu"); type = :ME)
        @named ball = MTK.FMIComponent(Val(2); fmu, type = :ME)
        @variables y(t) = 1.0
        # the import marks the states of an FMU with event indicators irreducible, so
        # simplification keeps them as unknowns and the events have somewhere to write to
        @mtkcompile sys = System([D(y) ~ -y], t; systems = [ball])
        @test ext.SII.variable_index(sys, ball.h) !== nothing
        @test ext.SII.variable_index(sys, ball.v) !== nothing
        prob = ODEProblem(sys, [ball.h => 1.0, ball.v => 0.0], (0.0, 3.0))
        sol = solve(prob, Tsit5(); abstol = 1.0e-8, reltol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)
        velocities = sol[ball.v]
        bounces = findall(i -> sol.t[i] == sol.t[i + 1], 1:(length(sol.t) - 1))
        @test length(bounces) >= 2
        @test all(i -> velocities[i] < 0 <= velocities[i + 1], bounces)

        # declaring a state an input of the compiled system takes it out of the unknowns
        # regardless, and a state event that cannot write its result back would silently do
        # nothing, so an FMU declaring event indicators still fails loudly at construction
        sys_input = MTK.mtkcompile(
            System([D(y) ~ -y], t; systems = [ball], name = :outer); inputs = [ball.v]
        )
        @test ext.SII.variable_index(sys_input, ball.v) === nothing
        par = only(
            p for p in parameters(sys_input) if MTK.hasmetadata(p, ext.FMUEventMetadata)
        )
        @test_throws "cannot be written back" MTK.getmetadata(
            par, CallbackConstructionHook
        )(sys_input, par)
    end
end

@testset "FMU event composition" begin
    # `save_positions = (true, true)` duplicates the time of every event
    event_indices(sol) = findall(i -> sol.t[i] == sol.t[i + 1], 1:(length(sol.t) - 1))
    # the events at which the FMU reversed this ball's velocity
    bounces(sol, velocities) = [
        i for i in event_indices(sol) if velocities[i] < 0 <= velocities[i + 1]
    ]

    @testset "FMU, native and user events in one solve, v$Ver, ME" for Ver in (2, 3)
        fmu = loadFMU(joinpath(FMU_DIR, "BouncingBall$Ver.fmu"); type = :ME)
        @named ball = MTK.FMIComponent(Val(Ver); fmu, type = :ME)
        @variables x(t) = 1.0
        @parameters rate = 1.0
        # an `ImperativeAffect` keeps the FMU out of the affect; see the equational-affect
        # testset below
        native = MTK.SymbolicContinuousCallback(
            [x ~ 0.5],
            MTK.ImperativeAffect(modified = (; rate)) do m, o, c, i
                return (; rate = 2 * m.rate)
            end
        )
        @mtkcompile sys = System(
            [D(x) ~ -rate * x], t; systems = [ball], continuous_events = [native]
        )
        user_times = Float64[]
        user_cb = SciMLBase.DiscreteCallback(
            (u, t, integrator) -> t >= 1.0,
            integrator -> (push!(user_times, integrator.t); nothing);
            save_positions = (false, false)
        )
        prob = ODEProblem(
            sys, [ball.h => 1.0, ball.v => 0.0], (0.0, 3.0); callback = user_cb
        )
        callbacks = prob.kwargs[:callback]
        # the FMU's state events and the native event, beside the FMU's finalization and
        # step callbacks and the user's
        @test length(callbacks.continuous_callbacks) == 2
        @test length(callbacks.discrete_callbacks) == 3
        @test count(cb -> cb === user_cb, callbacks.discrete_callbacks) == 1

        sol = solve(prob, Tsit5(); abstol = 1.0e-8, reltol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)
        # the FMU's events still bounce the ball
        @test minimum(sol[ball.h]) >= -1.0e-8
        @test length(bounces(sol, sol[ball.v])) >= 2
        # the native event fired: `x` decays from 1, crosses 0.5 at `log(2)` and decays
        # twice as fast from there
        @test sol.ps[rate] == 2.0
        @test sol[x][end] ≈ 0.5 * exp(-2 * (3.0 - log(2))) rtol = 1.0e-5
        # the user callback fired, and only where its own condition holds
        @test !isempty(user_times)
        @test minimum(user_times) >= 1.0
    end

    @testset "an equational native affect cannot differentiate the FMU without FMISensitivity, v2, ME" begin
        fmu = loadFMU(joinpath(FMU_DIR, "BouncingBall2.fmu"); type = :ME)
        @named ball = MTK.FMIComponent(Val(2); fmu, type = :ME)
        @variables x(t) = 1.0
        @discretes k(t) = 1.0
        native = MTK.SymbolicContinuousCallback(
            [x ~ 0.5] => [k ~ 2 * Pre(k)]; discrete_parameters = [k], iv = t
        )
        @mtkcompile sys = System(
            [D(x) ~ -k * x], t; systems = [ball], continuous_events = [native]
        )
        prob = ODEProblem(sys, [ball.h => 1.0, ball.v => 0.0], (0.0, 3.0))
        # an equational affect inherits the system's observed equations, which for a Model
        # Exchange FMU are the FMU call, and solves them implicitly. Differentiating that
        # solve differentiates the FMU, which is a black box without FMISensitivity.
        @test_throws "FMISensitivity" solve(
            prob, Tsit5(); abstol = 1.0e-8, reltol = 1.0e-8
        )
    end

    @testset "two FMU instances in one system, v$Ver, ME" for Ver in (2, 3)
        fmu = loadFMU(joinpath(FMU_DIR, "BouncingBall$Ver.fmu"); type = :ME)
        # one `FMIComponent` per instance: each builds its own wrapper, with its own buffers
        # and its own FMU instance
        @named ball1 = MTK.FMIComponent(Val(Ver); fmu, type = :ME)
        @named ball2 = MTK.FMIComponent(Val(Ver); fmu, type = :ME)
        @variables x(t) = 1.0
        @mtkcompile sys = System([D(x) ~ -x], t; systems = [ball1, ball2])
        prob = ODEProblem(
            sys,
            [ball1.h => 1.0, ball1.v => 0.0, ball2.h => 0.5, ball2.v => 0.0],
            (0.0, 3.0)
        )
        callbacks = prob.kwargs[:callback]
        # one state-event callback per FMU, and one finalization plus one step callback each
        @test length(callbacks.continuous_callbacks) == 2
        @test length(callbacks.discrete_callbacks) == 4

        sol = solve(prob, Tsit5(); abstol = 1.0e-8, reltol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)
        for (ball, h₀) in ((ball1, 1.0), (ball2, 0.5))
            @test minimum(sol[ball.h]) >= -1.0e-8
            bounced = bounces(sol, sol[ball.v])
            @test !isempty(bounced)
            # a ball dropped from `h₀` reaches the floor at `sqrt(2h₀/g)`, `g = 9.81`: a
            # callback reading the other instance's states would bounce both at once
            @test sol.t[first(bounced)] ≈ sqrt(2 * h₀ / 9.81) atol = 1.0e-6
        end
    end
end
