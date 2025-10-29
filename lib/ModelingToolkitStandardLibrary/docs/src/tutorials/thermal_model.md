# Heat Conduction Model

This example demonstrates the thermal response of two masses connected by a conducting element.
The two masses have the same heat capacity but different initial temperatures (`T1=100 [°C]`, `T2=0 [°C]`).
The mass with the higher temperature will cool off, while the mass with the lower temperature heats up.
They will each asymptotically approach the calculated temperature T_final_K that results
from dividing the total initial energy in the system by the sum of the heat capacities of each element.

```@example
using ModelingToolkitStandardLibrary.Thermal, ModelingToolkit, OrdinaryDiffEq, Plots
using ModelingToolkit: t_nounits as t

@mtkmodel HeatConductionModel begin
    @parameters begin
        C1 = 15
        C2 = 15
    end
    @components begin
        mass1 = HeatCapacitor(C = C1, T = 373.15)
        mass2 = HeatCapacitor(C = C2, T = 273.15)
        conduction = ThermalConductor(G = 10)
        Tsensor1 = TemperatureSensor()
        Tsensor2 = TemperatureSensor()
    end
    @equations begin
        connect(mass1.port, conduction.port_a)
        connect(conduction.port_b, mass2.port)
        connect(mass1.port, Tsensor1.port)
        connect(mass2.port, Tsensor2.port)
    end
end

@mtkcompile sys = HeatConductionModel()
prob = ODEProblem(sys, Pair[], (0, 5.0))
sol = solve(prob)

T_final_K = sol[(sys.mass1.T * sys.C1 + sys.mass2.T * sys.C2) / (sys.C1 + sys.C2)]

plot(title = "Thermal Conduction Demonstration")
plot!(sol, idxs = [sys.mass1.T, sys.mass2.T],
    labels = ["Mass 1 Temperature" "Mass 2 Temperature"])
plot!(sol.t, T_final_K, label = "Steady-State Temperature")
```
