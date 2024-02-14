include("electrical_components.jl")

@named source = ConstantVoltage(V = 10.0)
@named resistor = Resistor(R = 1.0)
@named inductor1 = Inductor(L = 1.0e-2)
@named inductor2 = Inductor(L = 2.0e-2)
@named ground = Ground()

eqs = [connect(source.p, resistor.p)
       connect(resistor.n, inductor1.p)
       connect(inductor1.n, inductor2.p)
       connect(source.n, inductor2.n)
       connect(inductor2.n, ground.g)]

@named ll_model = ODESystem(eqs, t)
ll_model = compose(ll_model, [source, resistor, inductor1, inductor2, ground])

@named source = ConstantVoltage(V = 10.0)
@named resistor1 = Resistor(R = 1.0)
@named resistor2 = Resistor(R = 1.0)
@named inductor1 = Inductor(L = 1.0e-2)
@named inductor2 = Inductor(L = 2.0e-2)
@named ground = Ground()

eqs = [connect(source.p, inductor1.p)
       connect(inductor1.n, resistor1.p)
       connect(inductor1.n, resistor2.p)
       connect(resistor1.n, resistor2.n)
       connect(resistor2.n, inductor2.p)
       connect(source.n, inductor2.n)
       connect(inductor2.n, ground.g)]
@named ll2_model = ODESystem(eqs, t)
ll2_model = compose(ll2_model, [source, resistor1, resistor2, inductor1, inductor2, ground])
