# Custom Component

In this tutorial, the creation of a custom component is demonstrated via the [Chua's circuit](https://en.wikipedia.org/wiki/Chua%27s_circuit).
The circuit is a simple circuit that shows chaotic behavior.
Except for a non-linear resistor, every other component already is part of `ModelingToolkitStandardLibrary.Electrical`.

First, we need to make some imports.

```@example components
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ModelingToolkitStandardLibrary.Electrical
using OrdinaryDiffEq
using Plots
```

## Custom Component

Now the custom component can be defined.
The Modelica implementation of the `NonlinearResistor` looks as follows:

```Modelica
model NonlinearResistor "Chua's resistor"
  extends Interfaces.OnePort;

  parameter SI.Conductance Ga "conductance in inner voltage range";
  parameter SI.Conductance Gb "conductance in outer voltage range";
  parameter SI.Voltage Ve "inner voltage range limit";
equation
  i = if (v < -Ve) then Gb*(v + Ve) - Ga*Ve else if (v > Ve) then Gb*(v - Ve) + Ga*Ve else Ga*v;
end NonlinearResistor;
```

this can almost be directly translated to the syntax of `ModelingToolkit`.

```@example components
@mtkmodel NonlinearResistor begin
    @extend OnePort()
    @parameters begin
        Ga
        Gb
        Ve
    end
    @equations begin
        i ~ ifelse(v < -Ve,
            Gb * (v + Ve) - Ga * Ve,
            ifelse(v > Ve,
                Gb * (v - Ve) + Ga * Ve,
                Ga * v))
    end
end
nothing # hide
```

### Explanation

Since the non-linear resistor is essentially a standard electrical component with two ports, we can extend from the `OnePort` component of the library.

```julia
@extend OnePort()
```

This extends `OnePort` and unpacks `v` and `i` variables.

It might be a good idea to create parameters for the constants of the `NonlinearResistor`.

```julia
@parameters begin
    Ga
    Gb
    Ve
end
```

This creates symbolic parameters with the name `Ga`, `Gb` and `Ve` whose default values are set from the function's arguments `Ga`, `Gb` and `Ve`, respectively.
This allows the user to `remake` the problem easily with different parameters or allow for auto-tuning or parameter optimization without having to do all the costly steps that may be involved with building and simplifying a model.
The non-linear (in this case piece-wise constant) equation for the current can be implemented using `ifelse`.

## Building the Model

The final model can now be created with the components from the library and the new custom component.

```@example components
@mtkmodel ChaoticAttractor begin
    @components begin
        inductor = Inductor(L = 18, i = 0)
        resistor = Resistor(R = 12.5e-3)
        conductor = Conductor(G = 0.565)
        capacitor1 = Capacitor(C = 10, v = 4)
        capacitor2 = Capacitor(C = 100, v = 0)
        non_linear_resistor = NonlinearResistor(
            Ga = -0.757576,
            Gb = -0.409091,
            Ve = 1
        )
        ground = Ground()
    end
    @equations begin
        connect(inductor.p, conductor.p)
        connect(conductor.n, non_linear_resistor.p)
        connect(capacitor1.p, conductor.n)
        connect(inductor.n, resistor.p)
        connect(conductor.p, capacitor2.p)
        connect(capacitor1.n, capacitor2.n, non_linear_resistor.n, resistor.n, ground.g)
    end
end
nothing # hide
```

## Simulating the Model

`@mtkcompile` builds a structurally simplified `ChaoticAttractor` model.
Since the initial voltage of the capacitors was already specified via `v` and the initial current of inductor via `i`, no initial condition is given and an empty pair is supplied.

```@example components
@mtkcompile sys = ChaoticAttractor()
prob = ODEProblem(sys, Pair[], (0, 5e4))
sol = solve(prob; saveat = 1.0)

plot(sol[sys.capacitor1.v], sol[sys.capacitor2.v], title = "Chaotic Attractor", label = "",
    ylabel = "C1 Voltage in V", xlabel = "C2 Voltage in V")
```

```@example components
plot(sol; idxs = [sys.capacitor1.v, sys.capacitor2.v, sys.inductor.i],
    labels = ["C1 Voltage in V" "C2 Voltage in V" "Inductor Current in A"])
```
