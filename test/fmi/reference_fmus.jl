using ModelingToolkit, FMI, FMIZoo, OrdinaryDiffEq, SciMLBase
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkit as MTK

# The FMI standard's Reference FMUs (https://github.com/modelica/Reference-FMUs), which
# FMIZoo downloads. They cover what the FMUs in `fmus/` do not: Model Exchange events,
# integer, boolean and enumeration variables, FMUs without states, and outputs that are
# states. Each result is compared with FMI.jl's own simulator. v0.0.39 is the newest
# release FMIZoo can download: later releases name their archive differently.
const REFERENCE_FMUS = "0.0.39"

function reference_fmu(model, version; type)
    return loadFMU(model, "ModelicaReferenceFMUs", REFERENCE_FMUS, "$version.0"; type)
end

const STEP = 0.01
const FMI3_INTERNAL_STEPS = Dict("Dahlquist" => 0.1, "Feedthrough" => 0.1, "Stair" => 0.2)

communication_step(model, version) = version == 3 ? get(FMI3_INTERNAL_STEPS, model, STEP) : STEP

function reference_component(model, version, type; name = :fmu)
    fmu = reference_fmu(model, version; type)
    kwargs = type == :CS ? (; communication_step_size = communication_step(model, version)) : (;)
    return MTK.FMIComponent(Val(version); fmu, type, name, fmu_initial_states = true, kwargs...)
end

# The values of the real variables `names` at each of `times`, from FMI.jl. Its
# Co-Simulation simulator records only the start and the end, so each time is its own run.
function fmi_values(model, version, type, names, times)
    return map(times) do time
        fmu = reference_fmu(model, version; type)
        kwargs = (; recordValues = names, saveat = [time], showProgress = false)
        sol = type == :ME ?
            simulateME(fmu, (0.0, time); solver = Tsit5(), kwargs...) :
            simulateCS(fmu, (0.0, time); dt = communication_step(model, version), kwargs...)
        @test sol.success
        collect(Float64, sol.values.saveval[end])
    end
end

# Model Exchange is integrated by different solvers on the two sides; Co-Simulation is
# stepped by the FMU itself, so it must agree exactly.
tolerance(type) = type == :ME ? (; rtol = 1.0e-6, atol = 1.0e-8) : (; rtol = 0.0, atol = 1.0e-12)

@testset "States and outputs match FMI.jl" begin
    # VanDerPol has two states; in its FMI 3 version the outputs are the states. BouncingBall bounces (a state event). Dahlquist checks Co-Simulation's last
    # step, at the end time. No initial conditions are given: the FMU sets them.
    cases = [
        ("VanDerPol", (2, 3), ["x0", "x1"]),
        ("BouncingBall", (2, 3), ["h", "v"]),
        ("Dahlquist", (2, 3), ["x"]),
    ]
    times = [0.5, 1.0]
    @testset "$model v$version $type" for (model, versions, names) in cases,
            version in versions, type in (:ME, :CS)

        @mtkcompile sys = reference_component(model, version, type)
        sol = solve(ODEProblem(sys, Dict(), (0.0, 1.0)), Tsit5())
        @test SciMLBase.successful_retcode(sol)
        expected = fmi_values(model, version, type, names, times)
        vars = [getproperty(sys, Symbol(n)) for n in names]
        for (i, time) in enumerate(times)
            @test [sol(time; idxs = v) for v in vars] ≈ expected[i] rtol = tolerance(type).rtol atol = tolerance(type).atol
        end
        if version == 3 && model == "VanDerPol"
            @test all(MTK.isoutput, vars)
        end
    end
end

@testset "Integer, boolean and enumeration ports" begin
    # Feedthrough copies each input to the output of the same type. It has no states.
    values = Dict(
        "Float64_continuous" => 1.5, "Int32" => 3.0, "Boolean" => 1.0, "Enumeration" => 2.0
    )
    @testset "v$version $type" for version in (2, 3), type in (:ME, :CS)
        @named fmu = reference_component("Feedthrough", version, type)
        @test MTK.isinput(fmu.Int32_input)
        @test MTK.isoutput(fmu.Boolean_output)
        inputs = Set(MTK.getname(v) for v in unknowns(fmu) if MTK.isinput(v))
        eqs = [
            getproperty(fmu, Symbol(n, "_input")) ~ get(values, n, 0.0)
                for n in (chopsuffix(string(i), "_input") for i in inputs)
        ]
        @mtkcompile sys = System(eqs, t; systems = [fmu])
        sol = solve(ODEProblem(sys, Dict(), (0.0, 1.0)), Tsit5())
        @test SciMLBase.successful_retcode(sol)
        for (n, value) in values
            @test sol[getproperty(sys.fmu, Symbol(n, "_output"))][end] == value
        end
    end
end

@testset "Time events and step-like outputs" begin
    # Stair counts up by one each second, from 1, through time events. It has no states.
    @testset "v$version $type" for version in (2, 3), type in (:ME, :CS)
        @mtkcompile sys = reference_component("Stair", version, type)
        offset = version == 3 && type == :CS ? 0.1 : 0.0
        prob = ODEProblem(sys, Dict(), (0.0, 3.5 + offset))
        sol = solve(prob, Tsit5())
        @test SciMLBase.successful_retcode(sol)
        @test [sol(time; idxs = sys.counter) for time in 0.5:1.0:3.5] == [1, 2, 3, 4]
        # Solving again, and over a shorter span, starts the count afresh.
        @test solve(prob, Tsit5())[sys.counter][end] == 4
        @test solve(remake(prob; tspan = (0.0, 1.5 + offset)), Tsit5())[sys.counter][end] == 2
    end
end

@testset "Co-Simulation over other time spans" begin
    @testset "v$version" for version in (2, 3)
        @mtkcompile sys = reference_component("VanDerPol", version, :CS)
        prob = ODEProblem(sys, Dict(), (0.0, 1.0))
        @test SciMLBase.successful_retcode(solve(prob, Tsit5()))
        # A longer span after a shorter one, and an end time between two steps.
        expected = fmi_values("VanDerPol", version, :CS, ["x0"], [2.0])
        longer = solve(remake(prob; tspan = (0.0, 2.0)), Tsit5())
        @test SciMLBase.successful_retcode(longer)
        @test longer[sys.x0][end] ≈ only(expected[end]) atol = 1.0e-12
        off_grid = remake(prob; tspan = (0.0, 1.005))
        if version == 3
            @test_throws r"FMI 3.*1\.005.*lastSuccessfulTime" solve(off_grid, Tsit5())
        else
            between = solve(off_grid, Tsit5())
            @test SciMLBase.successful_retcode(between)
            @test between.t[end] == 1.005
        end
    end
end

@testset "Initial states come from the FMU" begin
    @testset "Dahlquist $type" for type in (:ME, :CS)
        @mtkcompile sys = reference_component("Dahlquist", 2, type)
        @test solve(ODEProblem(sys, Dict(), (0.0, 1.0)), Tsit5())[sys.x][1] == 1.0
        # A value given to the problem still wins.
        given = solve(ODEProblem(sys, Dict(sys.x => 3.0), (0.0, 1.0)), Tsit5())
        @test given[sys.x][1] == 3.0
    end
    @testset "SimpleAdder $type" for type in (:ME, :CS)
        fmu = loadFMU(joinpath(@__DIR__, "fmus", "SimpleAdder.fmu"); type)
        kwargs = type == :CS ? (; communication_step_size = STEP) : (;)
        @named adder = MTK.FMIComponent(Val(2); fmu, type, fmu_initial_states = true, kwargs...)
        @mtkcompile sys = System([adder.a ~ 1.0, adder.b ~ 2.0], t; systems = [adder])
        sol = solve(ODEProblem(sys, Dict(), (0.0, 0.1)), Tsit5())
        @test SciMLBase.successful_retcode(sol)
        @test sol[sys.adder.c][1] == 1.0
    end
    @testset "a state with no start value" begin
        # SimpleArrayModel's states are calculated by the FMU when it is initialized.
        path = joinpath(@__DIR__, "fmus", "SimpleArrayModel.fmu")
        @named model = MTK.FMIComponent(Val(2); fmu = loadFMU(path; type = :ME), type = :ME, fmu_initial_states = true)
        sys = mtkcompile(model)
        sol = solve(ODEProblem(sys, Dict(), (0.0, 1.0)), Tsit5())
        @test SciMLBase.successful_retcode(sol)
        names = ["x[1,1]", "x[3,2]"]
        expected = simulateME(
            loadFMU(path; type = :ME), (0.0, 1.0);
            recordValues = names, saveat = [0.0, 1.0], showProgress = false
        ).values.saveval
        vars = [getproperty(sys, Symbol(n)) for n in names]
        @test [sol[v][1] for v in vars] ≈ collect(expected[1])
        @test [sol[v][end] for v in vars] ≈ collect(expected[end]) rtol = 1.0e-6
    end
end

@testset "Solving again after changing a parameter" begin
    @testset "$type" for type in (:ME, :CS)
        function component_system()
            fmu = loadFMU(joinpath(@__DIR__, "fmus", "SimpleAdder.fmu"); type)
            kwargs = type == :CS ? (; communication_step_size = STEP) : (;)
            @named adder = MTK.FMIComponent(Val(2); fmu, type, kwargs...)
            @variables z(t) = 0.0
            @mtkcompile sys = System(
                [adder.a ~ 1.0, adder.b ~ 2.0, D(z) ~ 1.0], t; systems = [adder]
            )
            return sys
        end
        sys = component_system()
        prob = ODEProblem(sys, Dict(sys.adder.c => 0.0), (0.0, 1.0))
        first_solve = solve(prob, Tsit5())
        @test SciMLBase.successful_retcode(first_solve)
        @test first_solve[sys.adder.out2][end] ≈ 8.0

        changed = solve(remake(prob; p = [sys.adder.value => 5.0]), Tsit5())
        fresh_sys = component_system()
        fresh = solve(
            ODEProblem(fresh_sys, Dict(fresh_sys.adder.c => 0.0, fresh_sys.adder.value => 5.0), (0.0, 1.0)),
            Tsit5()
        )
        @test SciMLBase.successful_retcode(changed)
        @test SciMLBase.successful_retcode(fresh)
        @test changed[sys.adder.out2][end] ≈ 16.0
        @test changed[sys.adder.out2][end] ≈ fresh[fresh_sys.adder.out2][end]
    end
end

@testset "Solvers that need derivatives" begin
    @testset "$model" for model in ("VanDerPol", "BouncingBall")
        @mtkcompile sys = reference_component(model, 3, :ME)
        prob = ODEProblem(sys, Dict(), (0.0, 1.0))
        stiff = solve(prob, Rodas5P(autodiff = AutoFiniteDiff()); reltol = 1.0e-8, abstol = 1.0e-10)
        @test SciMLBase.successful_retcode(stiff)
        @test stiff.stats.njacs > 0
        explicit = solve(prob, Tsit5(); reltol = 1.0e-8, abstol = 1.0e-10)
        @test stiff.u[end] ≈ explicit.u[end] rtol = 1.0e-6
    end
end
