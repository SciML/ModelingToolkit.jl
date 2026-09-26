# Tests of `linearize_hybrid` with components from DiscreteComponents and BlockComponents,
# simulated with SynchToolkit. These packages are registered in the DyadRegistry only, so this
# file is not part of the test groups run in CI. Run it in an environment that provides
# ModelingToolkit, DiscreteComponents, BlockComponents, SynchToolkit, ControlSystemsBase and
# OrdinaryDiffEqTsit5. The `DiscreteStateSpace` cases shift unscalarized array variables and
# were run with JuliaComputing/StateSelection.jl#161.
using ModelingToolkit, Test, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: unwrap
using DiscreteComponents
using DiscreteComponents: DiscretizationMethod
using BlockComponents.Continuous: FirstOrder
using BlockComponents.Math: Gain
using BlockComponents.Sources: Constant
using SynchToolkit
using OrdinaryDiffEqTsit5
import ControlSystemsBase as CS

# Convert a partition to a ControlSystemsBase state-space model.
to_ss(p) = p.Ts === nothing ? CS.ss(p.A, p.B, p.C, p.D) : CS.ss(p.A, p.B, p.C, p.D, p.Ts)

function freqresp_isapprox(G1, G2; w = exp10.(range(-2, 0.5, length = 20)), atol = 1.0e-8)
    return isapprox(CS.freqresp(G1, w), CS.freqresp(G2, w); atol)
end

const dt = 0.05

"""
Closed loop of a first-order plant and a discrete PI controller. With `with_reference = true`
the set point is a constant source, otherwise the set point input of the controller is left
unconnected so that it can serve as a linearization input.
"""
@component function PILoop(; name, with_reference = true, K = 2.0, Ti = 1.0)
    systems = @named begin
        plant = FirstOrder(k = 1, T = 1)
        sampler = Sampler()
        clock = DiscreteComponents.PeriodicClock(; dt)
        zoh = ZeroOrderHold()
        controller = DiscretePIDStandard(;
            K, Ti, Imethod = DiscretizationMethod.Forward(), with_D = false, y_max = Inf
        )
    end
    equations = [
        connect(controller.y, zoh.u)
        connect(zoh.y, plant.u)
        connect(plant.y, sampler.u)
        connect(sampler.y, controller.u_m, clock.y)
    ]
    if with_reference
        @named ref = Constant(k = 0.5)
        push!(systems, ref)
        push!(equations, connect(ref.y, controller.u_s))
    end
    return System(equations, t, [], []; name, systems)
end

@testset "Discrete PI controller loop" begin
    @named sim = PILoop()
    ssys = mtkcompile(sim; additional_passes = [SynchToolkit.compile_lustre])
    Tf = 3.0
    prob = ODEProblem(ssys, [ssys.plant.x => 0.0], (0.0, Tf))
    sol = solve(prob, Tsit5(); abstol = 1.0e-9, reltol = 1.0e-9)

    @named model = PILoop(; with_reference = false)
    model_nns = toggle_namespacing(model, false)
    hl = linearize_hybrid(
        model, [model_nns.controller.u_s], [model_nns.plant.y]; warn_missing_op = false
    )
    @test length(hl.partitions) == 2
    cont = hl.partitions[hl.continuous_index]
    disc = hl.partitions[2]
    @test disc.Ts == dt
    # The integrator state and the delayed integrator input of the forward Euler integrator.
    @test length(disc.unknowns) == 2
    @test cont.nu_user == 0
    @test cont.ny_user == 1
    @test disc.nu_user == 1
    @test disc.ny_user == 0
    @test length(hl.connections) == 2

    # Assemble the sampled-data loop from the set point to the plant output and compare its
    # step response with the simulation of the model.
    Pd = CS.c2d(to_ss(cont), dt)
    Cd = to_ss(disc)
    plant_to_ctrl = clock_boundary(hl, hl.continuous_index, 2)
    ctrl_to_plant = clock_boundary(hl, 2, hl.continuous_index)
    # The set point of the controller is the external input, the plant output the external
    # output.
    G = CS.feedback(
        Pd, Cd; Y1 = plant_to_ctrl.outputs, U2 = plant_to_ctrl.inputs,
        Y2 = ctrl_to_plant.outputs, U1 = ctrl_to_plant.inputs,
        W1 = 1:cont.nu_user, W2 = 1:disc.nu_user, Z1 = 1:cont.ny_user, pos_feedback = true
    )
    timevec = 0:dt:Tf
    res = CS.lsim(G, (x, t) -> [0.5], timevec)
    @test sol(timevec; idxs = ssys.plant.y).u ≈ res.y[:] rtol = 1.0e-5
end

@testset "DiscreteStateSpace" begin
    @testset "SISO" begin
        A = reshape([0.5], 1, 1)
        B = reshape([1.0], 1, 1)
        C = reshape([2.0], 1, 1)
        Dm = reshape([0.3], 1, 1)
        @component function SISOModel(; name)
            systems = @named begin
                sampler = Sampler()
                clock = DiscreteComponents.PeriodicClock(; dt)
                sys = DiscreteStateSpace(; nx = 1, nu = 1, ny = 1, A, B, C, D = Dm)
            end
            equations = [connect(sampler.y, sys.u[1], clock.y)]
            return System(equations, t, [], []; name, systems)
        end
        @named model = SISOModel()
        model_nns = toggle_namespacing(model, false)
        hl = linearize_hybrid(
            model, [model_nns.sampler.u], [model_nns.sys.y[1]]; warn_missing_op = false
        )
        disc = hl.partitions[2]
        @test disc.Ts == dt
        @test length(disc.unknowns) == 1
        @test freqresp_isapprox(to_ss(disc), CS.ss(A, B, C, Dm, dt))
        # The continuous partition only passes the input through to the sampler.
        cont = hl.partitions[1]
        @test isempty(cont.unknowns)
        @test cont.D == [1.0;;]
    end
    @testset "MIMO with array inputs and outputs" begin
        A = [0.5 0.1; -0.2 0.8]
        B = [1.0 0.0; 0.5 1.0]
        C = [1.0 0.0; 0.0 1.0]
        Dm = [0.0 0.1; 0.2 0.0]
        @component function MIMOModel(; name)
            systems = @named begin
                s1 = Sampler()
                s2 = Sampler()
                clock = DiscreteComponents.PeriodicClock(; dt)
                sys = DiscreteStateSpace(; nx = 2, nu = 2, ny = 2, A, B, C, D = Dm)
            end
            equations = [
                connect(s1.y, sys.u[1], clock.y)
                connect(s2.y, sys.u[2])
            ]
            return System(equations, t, [], []; name, systems)
        end
        @named model = MIMOModel()
        model_nns = toggle_namespacing(model, false)
        hl = linearize_hybrid(
            model, [model_nns.s1.u, model_nns.s2.u], [model_nns.sys.y]; warn_missing_op = false
        )
        disc = hl.partitions[2]
        @test length(disc.unknowns) == 2
        @test disc.ny_user == 2
        @test length(disc.inputs) == 2
        @test freqresp_isapprox(to_ss(disc), CS.ss(A, B, C, Dm, dt))
        @test sort(eigvals(disc.A)) ≈ sort(eigvals(A))
    end
end

@testset "MovingAverageFilter" begin
    N = 3
    @component function MAFModel(; name)
        systems = @named begin
            sampler = Sampler()
            clock = DiscreteComponents.PeriodicClock(; dt)
            filter = MovingAverageFilter(; N)
        end
        equations = [connect(sampler.y, filter.u, clock.y)]
        return System(equations, t, [], []; name, systems)
    end
    @named model = MAFModel()
    model_nns = toggle_namespacing(model, false)
    hl = linearize_hybrid(model, [model_nns.sampler.u], [model_nns.filter.y]; warn_missing_op = false)
    disc = hl.partitions[2]
    @test disc.Ts == dt
    @test length(disc.unknowns) == N - 1
    @test freqresp_isapprox(to_ss(disc), CS.tf([1.0, 1.0, 1.0] / N, [1.0, 0.0, 0.0], dt))
end

@testset "Unit delay whose output is consumed" begin
    # A one-sample delay `y ~ u(k - 1)` whose output enters another equation. The alias
    # elimination of the discrete partition then rebuilds the shifted form of the alias
    # target, which requires ModelingToolkitTearing to build it as a shift rather than as a
    # differential (fix of `var_derivative!` for discrete systems).
    systems = @named begin
        plant = FirstOrder(k = 1, T = 1)
        sampler = Sampler()
        clock = DiscreteComponents.PeriodicClock(; dt)
        delay = UnitDelay()
        gain = Gain(k = -0.5)
        zoh = ZeroOrderHold()
    end
    equations = [
        connect(plant.y, sampler.u)
        connect(sampler.y, delay.u, clock.y)
        connect(delay.y, gain.u)
        connect(gain.y, zoh.u)
        connect(zoh.y, plant.u)
    ]
    @named model = System(equations, t, [], []; systems)
    model_nns = toggle_namespacing(model, false)
    hl = linearize_hybrid(model, [], [model_nns.plant.y]; op = Dict(model_nns.plant.x => 0.0), warn_missing_op = false)
    disc = hl.partitions[2]
    @test length(disc.unknowns) == 1
    @test length(hl.connections) == 2
    # From the sampled plant output to the held controller output: a delayed gain.
    @test freqresp_isapprox(to_ss(disc), CS.tf(-0.5, [1, 0], dt))
end

@testset "Multi-rate loop" begin
    dt_fast = dt
    dt_slow = 4dt
    N = 4
    @component function MultiRateLoop(; name)
        systems = @named begin
            plant = FirstOrder(k = 1, T = 1)
            sampler = Sampler()
            fastclock = DiscreteComponents.PeriodicClock(; dt = dt_fast)
            filter = MovingAverageFilter(; N)
            lat = Latest()
            slowclock = DiscreteComponents.PeriodicClock(; dt = dt_slow)
            controller = DiscretePIDStandard(;
                K = 2.0, Ti = 1.0, Imethod = DiscretizationMethod.Forward(), with_D = false,
                y_max = Inf
            )
            zoh = ZeroOrderHold()
        end
        equations = [
            connect(plant.y, sampler.u)
            connect(sampler.y, filter.u, fastclock.y)
            connect(filter.y, lat.u)
            connect(lat.y, controller.u_m, slowclock.y)
            connect(controller.y, zoh.u)
            connect(zoh.y, plant.u)
        ]
        return System(equations, t, [], []; name, systems)
    end
    @named model = MultiRateLoop()
    model_nns = toggle_namespacing(model, false)
    hl = linearize_hybrid(
        model, [model_nns.controller.u_s], [model_nns.plant.y]; warn_missing_op = false
    )
    @test length(hl.partitions) == 3
    cont = hl.partitions[hl.continuous_index]
    fast = hl.partitions[findfirst(p -> p.Ts == dt_fast, hl.partitions)]
    slow = hl.partitions[findfirst(p -> p.Ts == dt_slow, hl.partitions)]
    @test length(hl.connections) == 3
    ifast = findfirst(p -> p.Ts == dt_fast, hl.partitions)
    islow = findfirst(p -> p.Ts == dt_slow, hl.partitions)
    @test length(clock_boundary(hl, hl.continuous_index, ifast).outputs) == 1
    @test length(clock_boundary(hl, ifast, islow).outputs) == 1
    @test length(clock_boundary(hl, islow, hl.continuous_index).outputs) == 1

    # The clock change is expressed with `Latest`, whose first argument is the fast signal.
    fast_to_slow = clock_boundary(hl, ifast, islow)
    term = slow.inputs[only(fast_to_slow.inputs)]
    @test operation(term) isa SynchToolkit.Latest
    @test isequal(fast.outputs[only(fast_to_slow.outputs)], unwrap(model_nns.lat.u))

    @test freqresp_isapprox(to_ss(fast), CS.tf([1.0, 1.0, 1.0, 1.0] / N, [1.0, 0.0, 0.0, 0.0], dt_fast))
    @test cont.A == [-1.0;;]
    @test length(slow.unknowns) == 2
    @test slow.nu_user == 1
end

"""
A plant with a cubic nonlinearity.
"""
@component function CubicPlant(; name)
    vars = @variables begin
        x(t) = 0.0
        u(t), [input = true]
        y(t), [output = true]
    end
    eqs = [D(x) ~ -x^3 + u, y ~ x]
    return System(eqs, t, vars, []; name)
end

@component function NonlinearLoop(; name, with_reference = true)
    systems = @named begin
        plant = CubicPlant()
        sampler = Sampler()
        clock = DiscreteComponents.PeriodicClock(; dt)
        zoh = ZeroOrderHold()
        controller = DiscretePIDStandard(;
            K = 2.0, Ti = 1.0, Imethod = DiscretizationMethod.Forward(), with_D = false,
            y_max = Inf
        )
    end
    equations = [
        connect(controller.y, zoh.u)
        connect(zoh.y, :plant_input, plant.u)
        connect(plant.y, sampler.u)
        connect(sampler.y, controller.u_m, clock.y)
    ]
    if with_reference
        @named ref = Constant(k = 1.0)
        push!(systems, ref)
        push!(equations, connect(ref.y, controller.u_s))
    end
    return System(equations, t, [], []; name, systems)
end

@testset "Operating point from a solution" begin
    @named sim = NonlinearLoop()
    ssys = mtkcompile(sim; additional_passes = [SynchToolkit.compile_lustre])
    prob = ODEProblem(ssys, [], (0.0, 4.0))
    sol = solve(prob, Tsit5(); abstol = 1.0e-9, reltol = 1.0e-9)
    t0 = 2.0
    @named model = NonlinearLoop(; with_reference = false)
    model_nns = toggle_namespacing(model, false)
    # The solution provides values of the observed variables as well, which must not enter
    # the initialization of the continuous partition as additional constraints.
    hl = @test_nowarn linearize_hybrid(
        model, [model_nns.controller.u_s], [model_nns.plant.y];
        op = LinearizationOpPoint(sol, t0; op = Dict(model_nns.controller.u_s => 1.0)),
        warn_missing_op = false
    )
    cont = hl.partitions[hl.continuous_index]
    x0 = sol(t0; idxs = ssys.plant.x)
    @test cont.x0 ≈ [x0]
    @test cont.A ≈ [-3x0^2;;]
    # The dictionary operating point gives the same result.
    hl2 = linearize_hybrid(
        model, [model_nns.controller.u_s], [model_nns.plant.y];
        op = Dict(model_nns.plant.x => x0, model_nns.controller.u_s => 1.0), warn_missing_op = false
    )
    @test hl2.partitions[hl2.continuous_index].A ≈ cont.A
end

@testset "Analysis points" begin
    @named model = NonlinearLoop()
    model_nns = toggle_namespacing(model, false)
    op = Dict(model_nns.plant.x => 0.5)
    hl = get_sensitivity(model, :plant_input; op, hybrid = true, warn_missing_op = false)
    @test hl isa HybridLinearization
    cont = hl.partitions[hl.continuous_index]
    @test cont.nu_user == 1
    @test cont.ny_user == 1
    @test cont.A ≈ [-3 * 0.5^2;;]
    hl = get_looptransfer(model, :plant_input; op, hybrid = true, warn_missing_op = false)
    @test hl isa HybridLinearization
    hl = linearize_hybrid(model, :plant_input, [model_nns.plant.y]; op, warn_missing_op = false)
    @test hl.partitions[hl.continuous_index].nu_user == 1
end
