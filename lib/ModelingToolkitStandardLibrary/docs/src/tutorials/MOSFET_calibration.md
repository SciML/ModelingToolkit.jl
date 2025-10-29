# MOSFET I-V Curves

In this example, first we'll demonstrate the I-V curves of the NMOS transistor model.
First of all, we construct a circuit using the NMOS transistor. We'll need to import ModelingToolkit and the Electrical standard library that holds the transistor models.

```@example NMOS
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Blocks: Constant
using OrdinaryDiffEq
using Plots
```

Here we just connect the source pin to ground, the drain pin to a voltage source named `Vcc`, and the gate pin to a voltage source named `Vb`.

```@example NMOS
@mtkmodel SimpleNMOSCircuit begin
    @components begin
        Q1 = NMOS()
        Vcc = Voltage()
        Vb = Voltage()
        ground = Ground()

        Vcc_const = Constant(k = V_cc)
        Vb_const = Constant(k = V_b)
    end

    @parameters begin
        V_cc = 5.0
        V_b = 3.5
    end
    @equations begin
        #voltage sources
        connect(Vcc_const.output, Vcc.V)
        connect(Vb_const.output, Vb.V)

        #ground connections
        connect(Vcc.n, Vb.n, ground.g, Q1.s)

        #other stuff
        connect(Vcc.p, Q1.d)
        connect(Vb.p, Q1.g)
    end
end

@mtkbuild sys = SimpleNMOSCircuit(V_cc = 5.0, V_b = 3.5)

prob = ODEProblem(sys, Pair[], (0.0, 10.0))
sol = solve(prob)
```

Now to make sure that the transistor model is working like it's supposed to, we can examine the plots of the drain-source voltage vs. the drain current, otherwise knowns as the I-V curve of the transistor.

```@example NMOS
v_cc_list = collect(0.05:0.1:10.0)

I_D_list = []
I_D_lists = []

for V_b in [1.0, 1.4, 1.8, 2.2, 2.6]
    I_D_list = []
    for V_cc in v_cc_list
        @mtkbuild sys = SimpleNMOSCircuit(V_cc = V_cc, V_b = V_b)
        prob = ODEProblem(sys, Pair[], (0.0, 10.0))
        sol = solve(prob)
        push!(I_D_list, sol[sys.Q1.d.i][1])
    end
    push!(I_D_lists, I_D_list)
end

reduce(hcat, I_D_lists)
plot(v_cc_list, I_D_lists, title = "NMOS IV Curves",
    label = ["V_GS: 1.0 V" "V_GS: 1.4 V" "V_GS: 1.8 V" "V_GS: 2.2 V" "V_GS: 2.6 V"],
    xlabel = "Drain-Source Voltage (V)", ylabel = "Drain Current (A)")
```

We can see that we get exactly what we would expect: as the drain-source voltage increases, the drain current increases, until the the transistor gets in to the saturation region of operation.
Then the only increase in drain current is due to the channel-length modulation effect. Additionally, we can see that the maximum current reached increases as the gate voltage increases.
