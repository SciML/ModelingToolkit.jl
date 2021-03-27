include("electrical_components.jl")

@named source = ConstantVoltage(V=10.0)
@named resistor = Resistor(R=1.0)
@named inductor1 = Inductor(L=8.0e-9)
@named inductor2 = Inductor(L=2.0e-9)
@named ground = Ground()

eqs = [
       connect_pins(source.p, resistor.p)
       connect_pins(resistor.n, inductor1.p)
       connect_pins(inductor1.n, inductor2.p)
       connect_pins(source.n, inductor2.n, ground.g)
      ]

@named ll_model = ODESystem(eqs, t, systems=[source, resistor, inductor1, inductor2, ground])
