using ModelingToolkitStandardLibrary.Thermal, ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitStandardLibrary.Blocks: Constant, Step
using OrdinaryDiffEq: ReturnCode.Success

# Test HeatCapacitor, TemperatureSensor, RelativeTemperatureSensor, FixedTemperature
@testset "Heat systems" begin
    T, C, G = 10.0, 10.0, 10.0
    @variables final_T(t)
    @named mass1 = HeatCapacitor(C = C)
    @named mass2 = HeatCapacitor(C = C)
    @named th_conductor = ThermalConductor(G = G)
    @named reltem_sensor = RelativeTemperatureSensor()
    @named T_sensor1 = TemperatureSensor()
    @named T_sensor2 = TemperatureSensor()
    @named tem_src = FixedTemperature(T = T)

    @info "Building a single-body system..."
    eqs = [connect(mass1.port, th_conductor.port_a)
           connect(th_conductor.port_b, reltem_sensor.port_a)
           connect(reltem_sensor.port_b, tem_src.port)]
    @named h1 = System(eqs, t, systems = [mass1, reltem_sensor, tem_src, th_conductor])
    sys = mtkcompile(h1)

    u0 = [mass1.T => 2]
    prob = ODEProblem(sys, u0, (0, 2.0))
    sol = solve(prob, Tsit5())

    # Check if Relative temperature sensor reads the temperature of heat capacitor
    # when connected to a thermal conductor and a fixed temperature source
    @test SciMLBase.successful_retcode(sol)
    @test sol[reltem_sensor.T.u] + sol[tem_src.port.T] ==
          sol[mass1.T] + sol[th_conductor.dT]

    @info "Building a two-body system..."
    eqs = [connect(T_sensor1.port, mass1.port, th_conductor.port_a)
           connect(th_conductor.port_b, mass2.port, T_sensor2.port)
           final_T ~ (mass1.C * mass1.T + mass2.C * mass2.T) /
                     (mass1.C + mass2.C)]
    @named h2 = System(eqs, t, [final_T], [],
        systems = [mass1, mass2, T_sensor1, T_sensor2, th_conductor])
    sys = mtkcompile(h2)

    u0 = [mass1.T => 1.0
          mass2.T => 10.0]
    prob = ODEProblem(sys, u0, (0, 3.0))
    sol = solve(prob, Tsit5())

    @test SciMLBase.successful_retcode(sol)
    m1, m2 = sol.u[end]
    @test m1≈m2 atol=1e-1
    @test sol[T_sensor1.T.u] == sol[sys.mass1.T]
    @test sol[T_sensor2.T.u] == sol[sys.mass2.T]
end

# Test HeatFlowSensor, FixedHeatFlow, ThermalResistor, ThermalConductor
@testset "Heat flow system" begin
    C, G, R = 10, 10, 10
    @named flow_src = FixedHeatFlow(Q_flow = 50, alpha = 100)
    @named mass1 = HeatCapacitor(C = C)
    @named hf_sensor1 = HeatFlowSensor()
    @named hf_sensor2 = HeatFlowSensor()
    @named th_conductor = ThermalConductor(G = G)
    @named th_resistor = ThermalResistor(R = R)
    @named th_ground = FixedTemperature(T = 0)

    @info "Building a heat-flow system..."
    eqs = [connect(mass1.port, th_resistor.port_a, th_conductor.port_a)
           connect(th_conductor.port_b, flow_src.port, hf_sensor1.port_a,
               hf_sensor2.port_a)
           connect(th_resistor.port_b, hf_sensor1.port_b, hf_sensor2.port_b,
               th_ground.port)]
    @named h2 = System(eqs, t,
        systems = [mass1, hf_sensor1, hf_sensor2,
            th_resistor, flow_src, th_ground, th_conductor])
    sys = mtkcompile(h2)

    u0 = [mass1.T => 10.0]
    prob = ODEProblem(sys, u0, (0, 3.0))
    sol = solve(prob, Tsit5())

    @test SciMLBase.successful_retcode(sol)
    @test sol[th_conductor.dT] .* G == sol[th_conductor.Q_flow]
    @test sol[th_conductor.Q_flow] ≈ sol[hf_sensor1.Q_flow.u] + sol[flow_src.port.Q_flow]

    @test sol[mass1.T] == sol[th_resistor.port_a.T]
    @test sol[th_resistor.dT] ./ R ≈ sol[th_resistor.Q_flow]
end

# Test ConvectiveConductor, BodyRadiation
@testset "Radiator system" begin
    T_gas, T_coolant = 1000, 10
    R_wall = 10
    G = 0.04
    σ = 5.6703744191844294e-8 # Stefan-Boltzmann constant

    @named base = ThermalResistor(R = R_wall)
    @named gas_tem = FixedTemperature(T = T_gas)
    @named coolant_tem = FixedTemperature(T = T_coolant)
    @named radiator = BodyRadiation(G = G)
    @named dissipator = ConvectiveConductor(G = 10)
    @named mass = HeatCapacitor(C = 10)

    @info "Building a radiator..."
    eqs = [connect(gas_tem.port, radiator.port_a, base.port_a, dissipator.solid, mass.port)
           connect(coolant_tem.port, base.port_b, radiator.port_b, dissipator.fluid)]
    @named rad = System(eqs, t,
        systems = [
            base,
            gas_tem,
            radiator,
            dissipator,
            coolant_tem,
            mass
        ])
    sys = mtkcompile(rad)

    u0 = [mass.T => T_gas]
    prob = ODEProblem(sys, u0, (0, 3.0))
    sol = solve(prob, Rodas4())

    @test SciMLBase.successful_retcode(sol)
    @test sol[dissipator.dT] == sol[radiator.port_a.T] - sol[radiator.port_b.T]
    rad_Q_flow = G * σ * (T_gas^4 - T_coolant^4)
    @test sol[radiator.Q_flow] == fill(rad_Q_flow, length(sol[radiator.Q_flow]))
end

@testset "Thermal Collector" begin
    @named flow_src = FixedHeatFlow(Q_flow = 50, alpha = 100)
    @named hf_sensor = HeatFlowSensor()
    @named collector = ThermalCollector(m = 2)
    @named th_resistor = ThermalResistor(R = 10)
    @named tem_src = FixedTemperature(T = 10)
    @named mass = HeatCapacitor(C = 10)

    @info "Building a heat collector..."
    eqs = [connect(flow_src.port, collector.port_a1, th_resistor.port_a)
           connect(tem_src.port, collector.port_a2)
           connect(hf_sensor.port_a, collector.port_b)
           connect(hf_sensor.port_b, mass.port, th_resistor.port_b)]
    @named coll = System(eqs, t,
        systems = [hf_sensor, flow_src, tem_src,
            collector, th_resistor, mass])
    sys = mtkcompile(coll)

    prob = ODEProblem(sys, [], (0, 3.0))
    sol = solve(prob, Rodas4())

    @test SciMLBase.successful_retcode(sol)
    @test sol[collector.port_b.Q_flow] + sol[collector.port_a1.Q_flow] +
          sol[collector.port_a2.Q_flow] ==
          zeros(length(sol[collector.port_b.Q_flow]))
    @test sol[collector.port_b.T] == sol[collector.port_a1.T] == sol[collector.port_a2.T]
end

@testset "FixedHeatFlow with alpha=0.0 test" begin
    @mtkmodel TestModel begin
        @components begin
            temp = FixedTemperature(T = 300)
            heatflow = FixedHeatFlow(Q_flow = -1.0)
            wall = ThermalResistor(R = 1)
        end

        @equations begin
            connect(temp.port, wall.port_a)
            connect(wall.port_b, heatflow.port)
        end
    end

    @info "Building a FixedHeatFlow with alpha=0.0"
    @mtkcompile test_model = TestModel()
    prob = ODEProblem(test_model, Pair[], (0, 10.0))
    sol = solve(prob)

    heat_flow = sol[test_model.heatflow.port.Q_flow]

    @test SciMLBase.successful_retcode(sol) # Ensure the simulation is successful
    @test all(isapprox.(heat_flow, 1.0, rtol = 1e-6)) # Heat flow value should be equal to the fixed value defined
end
