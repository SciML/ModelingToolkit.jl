include("electrical_components.jl")

R = 1.0
C = 1.0
V = 1.0
@named resistor = Resistor(R = R)
@named capacitor = Capacitor(C = C)
@named source = ConstantVoltage(V = V)
@named ground = Ground()

rc_eqs = [connect(source.p, resistor.p)
          connect(resistor.n, capacitor.p)
          connect(capacitor.n, source.n)
          connect(capacitor.n, ground.g)]

@named rc_model = ODESystem(rc_eqs, t)
rc_model = compose(rc_model, [resistor, capacitor, source, ground])
