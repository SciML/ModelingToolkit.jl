using ModelingToolkitStandardLibrary.Thermal, ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq: ReturnCode.Success

# Modelica example
@testset "demo" begin
    @named mass1 = HeatCapacitor(C = 15, T = 373.15)
    @named mass2 = HeatCapacitor(C = 15, T = 273.15)
    @named conduction = ThermalConductor(G = 10)
    @named Tsensor1 = TemperatureSensor()
    @named Tsensor2 = TemperatureSensor()

    connections = [
        connect(mass1.port, conduction.port_a),
        connect(conduction.port_b, mass2.port),
        connect(mass1.port, Tsensor1.port),
        connect(mass2.port, Tsensor2.port)
    ]

    @named model = System(connections, t,
        systems = [mass1, mass2, conduction, Tsensor1, Tsensor2])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, [], (0, 3.0))
    sol = solve(prob, Tsit5())
    @test SciMLBase.successful_retcode(sol)
end
