using ModelingToolkit, Test
using ModelingToolkitStandardLibrary.Thermal
using ModelingToolkitStandardLibrary.Blocks

# https://doc.modelica.org/Modelica%204.0.0/Resources/helpWSM/Modelica/Modelica.Thermal.HeatTransfer.Examples.Motor.html

@testset "Thermal Motor Demo" begin
    k2c(T) = T - 273.15

    @mtkmodel ThermalMotor begin
        @parameters begin
            T_amb = 293.15
        end
        @components begin
            windingLosses = PrescribedHeatFlow(T_ref = k2c(95), alpha = 3.03e-3)
            winding = HeatCapacitor(C = 2500, T = T_amb)
            T_winding = TemperatureSensor()
            winding2core = ThermalConductor(G = 10)
            coreLosses = PrescribedHeatFlow()
            core = HeatCapacitor(C = 25000, T = T_amb)
            T_core = TemperatureSensor()
            convection = ConvectiveConductor(G = 25)
            environment = PrescribedTemperature()
            amb = Constant(k = T_amb)
            core_losses_const = Constant(k = 500)
            winding_losses = Step(height = 900, offset = 100, start_time = 360,
                duration = Inf, smooth = false)
        end
        @equations begin
            connect(windingLosses.port, winding.port)
            connect(coreLosses.port, core.port)
            connect(winding.port, winding2core.port_a)
            connect(winding2core.port_b, core.port)
            connect(winding.port, T_winding.port)
            connect(core.port, T_core.port)
            connect(winding2core.port_b, convection.solid)
            connect(convection.fluid, environment.port)
            connect(amb.output, environment.T)
            connect(winding_losses.output, windingLosses.Q_flow)
            connect(core_losses_const.output, coreLosses.Q_flow)
        end
    end

    @mtkcompile motor = ThermalMotor()
    prob = ODEProblem(motor, Pair[], (0, 720.0))
    sol = solve(prob)

    # plot(sol; vars=[T_winding.T, T_core.T])
    @test SciMLBase.successful_retcode(sol)
    @test sol[motor.T_winding.T.u] == sol[motor.winding.T]
    @test sol[motor.T_core.T.u] == sol[motor.core.T]
    @test sol[-motor.core.port.Q_flow] ≈
          sol[motor.coreLosses.port.Q_flow + motor.convection.solid.Q_flow + motor.winding2core.port_b.Q_flow]
    @test sol[motor.T_winding.T.u][end] >= 500 # not good but better than nothing
    @test sol[motor.T_core.T.u] <= sol[motor.T_winding.T.u]
end
