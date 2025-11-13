using ModelingToolkit, Test
using ModelingToolkitStandardLibrary.Thermal
using ModelingToolkitStandardLibrary.Blocks

# Tests ConvectiveResistor and includes FixedTemperature and ThermalResistor

@testset "Piston cylinder wall" begin
    @info "Building a piston-cylinder..."
    @mtkmodel Piston begin
        @parameters begin
            # ᵧ -> gas and ᵪ -> coolant
            Tᵧ = 1000, [description = "Temperature of gas"]
            Tᵪ = 10, [description = "Temperature of coolant"]
            # R = 1/h; h is convection co-efficient
            Rᵧ = 50e-4, [description = "Thermal resistance of gas"]
            Rᵪ = 10e-4, [description = "Thermal resistance of coolant"]
            R_wall = 1.5e-4
        end
        @components begin
            coolant = ConvectiveResistor(R = Rᵪ)
            gas = ConvectiveResistor(R = Rᵧ)
            wall = ThermalResistor(R = R_wall)
            gas_tem = FixedTemperature(T = Tᵧ)
            coolant_tem = FixedTemperature(T = Tᵪ)
        end
        @equations begin
            connect(gas_tem.port, gas.solid)
            connect(gas.fluid, wall.port_a)
            connect(wall.port_b, coolant.fluid)
            connect(coolant.solid, coolant_tem.port)
        end
    end

    @mtkcompile piston = Piston()

    prob = ODEProblem(piston, [], (0, 3.0))
    sol = solve(prob)

    # Heat-flow-rate is equal in magnitude
    # and opposite in direction
    @test SciMLBase.successful_retcode(sol)
    # The initial value doesn't add up to absolute zero, while the rest do. To avoid
    # tolerance on the latter, the test is split in two parts.
    @test sol[piston.gas.Q_flow][1]≈-sol[piston.coolant.Q_flow][1] rtol=1e-6
    @test sol[piston.gas.Q_flow][2:end]≈-sol[piston.coolant.Q_flow][2:end] rtol=1e-6
end
