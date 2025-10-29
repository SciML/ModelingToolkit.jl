using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Thermal
import ModelingToolkitStandardLibrary
using ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq: ReturnCode.Success
# using Plots

@testset "DC motor" begin
    f = 0.01
    k = 0.5
    R = 0.5
    tau_L_step = -3
    V_step = 10

    @mtkmodel DCMotor begin
        @structural_parameters begin
            R = 0.5
            L = 4.5e-3
            k = 0.5
            J = 0.02
            f = 0.01
            V_step = 10
            tau_L_step = -3
        end
        @components begin
            ground = Ground()
            source = Voltage()
            voltage_step = Blocks.Step(height = V_step, start_time = 0)
            R1 = Resistor(R = R)
            L1 = Inductor(L = L, i = 0.0)
            emf = EMF(k = k)
            fixed = Fixed()
            load = Torque()
            load_step = Blocks.Step(height = tau_L_step, start_time = 3)
            inertia = Inertia(J = J)
            friction = Damper(d = f)
        end
        @equations begin
            connect(fixed.flange, emf.support, friction.flange_b)
            connect(emf.flange, friction.flange_a, inertia.flange_a)
            connect(inertia.flange_b, load.flange)
            connect(load_step.output, load.tau)
            connect(voltage_step.output, source.V)
            connect(source.p, R1.p)
            connect(R1.n, L1.p)
            connect(L1.n, emf.p)
            connect(emf.n, source.n, ground.g)
        end
    end

    @mtkcompile dc_motor = DCMotor(; f, k, R, V_step)

    prob = ODEProblem(dc_motor, unknowns(dc_motor) .=> 0.0, (0, 6.0))
    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    # EMF equations
    @test -0.5 .* sol[dc_motor.emf.i] == sol[dc_motor.emf.flange.tau]
    @test sol[dc_motor.emf.v] == 0.5 .* sol[dc_motor.emf.w]
    # test steady-state values
    dc_gain = [f/(k^2 + f * R) k/(k^2 + f * R); k/(k^2 + f * R) -R/(k^2 + f * R)]
    idx_t = findfirst(sol.t .> 2.5)
    @test sol[dc_motor.inertia.w][idx_t]≈(dc_gain * [V_step; 0])[2] rtol=1e-3
    @test sol[dc_motor.emf.i][idx_t]≈(dc_gain * [V_step; 0])[1] rtol=1e-3
    idx_t = findfirst(sol.t .> 5.5)
    @test sol[dc_motor.inertia.w][idx_t]≈(dc_gain * [V_step; -tau_L_step])[2] rtol=1e-3
    @test sol[dc_motor.emf.i][idx_t]≈(dc_gain * [V_step; -tau_L_step])[1] rtol=1e-3

    @test_skip begin
        prob = DAEProblem(dc_motor, D.(unknowns(dc_motor)) .=> 0.0, (0, 6.0))
        sol = solve(prob, DFBDF())
        @test sol.retcode == Success
        # EMF equations
        @test -0.5 .* sol[dc_motor.emf.i] == sol[dc_motor.emf.flange.tau]
        @test sol[dc_motor.emf.v] == 0.5 .* sol[dc_motor.emf.w]
        # test steady-state values
        dc_gain = [f/(k^2 + f * R) k/(k^2 + f * R); k/(k^2 + f * R) -R/(k^2 + f * R)]
        idx_t = findfirst(sol.t .> 2.5)
        @test sol[dc_motor.inertia.w][idx_t]≈(dc_gain * [V_step; 0])[2] rtol=1e-3
        @test sol[dc_motor.emf.i][idx_t]≈(dc_gain * [V_step; 0])[1] rtol=1e-3
        idx_t = findfirst(sol.t .> 5.5)
        @test sol[dc_motor.inertia.w][idx_t]≈(dc_gain * [V_step; -tau_L_step])[2] rtol=1e-3
        @test sol[dc_motor.emf.i][idx_t]≈(dc_gain * [V_step; -tau_L_step])[1] rtol=1e-3
    end
    # p1 = Plots.plot(sol, vars=[inertia.w], ylabel="Angular Vel. in rad/s", label="")
    # p2 = Plots.plot(sol, vars=[emf.i], ylabel="Current in A", label="")
    # Plots.plot(p1, p2, layout=(2,1))
    # Plots.savefig("dc_motor.png")
end

@testset "DC motor with speed sensor" begin
    R = 0.5
    k = 0.5
    f = 0.01
    V_step = 10
    tau_L_step = -3

    @mtkmodel DCMotorWithSpeedSensor begin
        @structural_parameters begin
            R = 0.5
            L = 4.5e-3
            k = 0.5
            J = 0.02
            f = 0.01
            V_step = 10
            tau_L_step = -3
        end
        @components begin
            ground = Ground()
            source = Voltage()
            voltage_step = Blocks.Step(height = V_step, start_time = 0)
            R1 = Resistor(R = R)
            L1 = Inductor(L = L, i = 0.0)
            emf = EMF(k = k)
            fixed = Fixed()
            load = Torque()
            load_step = Blocks.Step(height = tau_L_step, start_time = 3)
            inertia = Inertia(J = J)
            friction = Damper(d = f)
            speed_sensor = SpeedSensor()
        end
        @equations begin
            connect(fixed.flange, emf.support, friction.flange_b)
            connect(emf.flange, friction.flange_a, inertia.flange_a)
            connect(inertia.flange_b, load.flange)
            connect(inertia.flange_b, speed_sensor.flange)
            connect(load_step.output, load.tau)
            connect(voltage_step.output, source.V)
            connect(source.p, R1.p)
            connect(R1.n, L1.p)
            connect(L1.n, emf.p)
            connect(emf.n, source.n, ground.g)
        end
    end

    @mtkcompile sys = DCMotorWithSpeedSensor(; f, k, R, V_step, tau_L_step)

    prob = ODEProblem(sys, unknowns(sys) .=> 0.0, (0, 6.0))
    sol = solve(prob, Rodas4())

    @test sol.retcode == Success
    # EMF equations
    @test -0.5 .* sol[sys.emf.i] == sol[sys.emf.flange.tau]
    @test sol[sys.emf.v] == 0.5 .* sol[sys.emf.w]

    # test steady-state values
    dc_gain = [f/(k^2 + f * R) k/(k^2 + f * R); k/(k^2 + f * R) -R/(k^2 + f * R)]
    idx_t = findfirst(sol.t .> 2.5)
    @test sol[sys.inertia.w][idx_t]≈(dc_gain * [V_step; 0])[2] rtol=1e-3
    @test sol[sys.emf.i][idx_t]≈(dc_gain * [V_step; 0])[1] rtol=1e-3
    idx_t = findfirst(sol.t .> 5.5)
    @test sol[sys.inertia.w][idx_t]≈(dc_gain * [V_step; -tau_L_step])[2] rtol=1e-3
    @test sol[sys.emf.i][idx_t]≈(dc_gain * [V_step; -tau_L_step])[1] rtol=1e-3
    @test all(sol[sys.inertia.w] .== sol[sys.speed_sensor.w.u])

    @test_skip begin
        prob = DAEProblem(sys, D.(unknowns(sys)) .=> 0.0, (0, 6.0))
        sol = solve(prob, DFBDF())
        @test sol.retcode == Success
        # EMF equations
        @test -0.5 .* sol[sys.emf.i] == sol[sys.emf.flange.tau]
        @test sol[sys.emf.v] == 0.5 .* sol[sys.emf.w]
        # test steady-state values
        dc_gain = [f/(k^2 + f * R) k/(k^2 + f * R); k/(k^2 + f * R) -R/(k^2 + f * R)]
        idx_t = findfirst(sol.t .> 2.5)
        @test sol[sys.inertia.w][idx_t]≈(dc_gain * [V_step; 0])[2] rtol=1e-3
        @test sol[sys.emf.i][idx_t]≈(dc_gain * [V_step; 0])[1] rtol=1e-3
        idx_t = findfirst(sol.t .> 5.5)
        @test sol[sys.inertia.w][idx_t]≈(dc_gain * [V_step; -tau_L_step])[2] rtol=1e-3
        @test sol[sys.emf.i][idx_t]≈(dc_gain * [V_step; -tau_L_step])[1] rtol=1e-3
        #
        @test all(sol[sys.inertia.w] .== sol[sys.speed_sensor.w.u])
    end
end

@testset "Electrical Heating Circuit" begin
    @mtkmodel ElHeatingCircuit begin
        @components begin
            ground = Ground()
            source = Voltage()
            voltage_sine = Blocks.Sine(amplitude = 220, frequency = 1)
            heating_resistor = Resistor(R = 100, alpha = 1e-3,
                T_ref = 293.15, T_dep = true)
            thermal_conductor = ThermalConductor(G = 50)
            env = FixedTemperature(T = 273.15 + 20)
        end
        @equations begin
            connect(source.n, ground.g, heating_resistor.n)
            connect(source.p, heating_resistor.p)
            connect(voltage_sine.output, source.V)
            connect(heating_resistor.heat_port, thermal_conductor.port_a)
            connect(thermal_conductor.port_b, env.port)
        end
    end

    @mtkcompile sys = ElHeatingCircuit()

    prob = ODEProblem(sys, [], (0, 6.0); guesses = [sys.heating_resistor.i => 0.0])
    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test sol[sys.source.v * sys.source.i] == -sol[sys.env.port.Q_flow]
end
