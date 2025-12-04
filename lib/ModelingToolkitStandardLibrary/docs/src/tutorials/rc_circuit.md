# RC Circuit Model

This tutorial is a simplified version of the [RC circuit tutorial in the
`ModelingToolkit.jl` documentation](https://docs.sciml.ai/ModelingToolkit/stable/tutorials/acausal_components/).
In that tutorial, the full RC circuit is built from scratch. Here, we will use the
components of the `Electrical` model in the ModelingToolkit Standard Library to simply
connect pre-made components and simulate the model.

```@example
using ModelingToolkit, OrdinaryDiffEq, Plots
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Blocks: Constant
using ModelingToolkit: t_nounits as t

@mtkmodel RC begin
    @parameters begin
        R = 1.0
        C = 1.0
        V = 1.0
    end
    @components begin
        resistor = Resistor(R = R)
        capacitor = Capacitor(C = C, v = 0.0)
        source = Voltage()
        constant = Constant(k = V)
        ground = Ground()
    end
    @equations begin
        connect(constant.output, source.V)
        connect(source.p, resistor.p)
        connect(resistor.n, capacitor.p)
        connect(capacitor.n, source.n, ground.g)
    end
end

@mtkcompile sys = RC()
prob = ODEProblem(sys, Pair[], (0, 10.0))
sol = solve(prob)

plot(sol, idxs = [sys.capacitor.v, sys.resistor.i],
    title = "RC Circuit Demonstration",
    labels = ["Capacitor Voltage" "Resistor Current"])
```
