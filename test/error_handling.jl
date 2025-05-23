using Test
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkit: ExtraVariablesSystemException, ExtraEquationsSystemException

using ModelingToolkitStandardLibrary.Electrical

function UnderdefinedConstantVoltage(; name, V = 1.0)
    val = V
    @named p = Pin()
    @named n = Pin()
    @parameters V
    eqs = [
        V ~ p.v - n.v        # Remove equation
    # 0 ~ p.i + n.i
    ]
    System(eqs, t, [], [V], systems = [p, n], defaults = Dict(V => val), name = name)
end

function OverdefinedConstantVoltage(; name, V = 1.0, I = 1.0)
    val = V
    val2 = I
    @named p = Pin()
    @named n = Pin()
    @parameters V I
    eqs = [V ~ p.v - n.v
           # Overdefine p.i and n.i
           n.i ~ I
           p.i ~ I]
    System(eqs, t, [], [V, I], systems = [p, n], defaults = Dict(V => val, I => val2),
        name = name)
end

R = 1.0
C = 1.0
V = 1.0
@named resistor = Resistor(R = R)
@named capacitor = Capacitor(C = C)
@named source = UnderdefinedConstantVoltage(V = V)

rc_eqs = [connect(source.p, resistor.p)
          connect(resistor.n, capacitor.p)
          connect(capacitor.n, source.n)]

@named rc_model = System(rc_eqs, t, systems = [resistor, capacitor, source])
@test_throws ModelingToolkit.ExtraVariablesSystemException mtkcompile(rc_model)

@named source2 = OverdefinedConstantVoltage(V = V, I = V / R)
rc_eqs2 = [connect(source2.p, resistor.p)
           connect(resistor.n, capacitor.p)
           connect(capacitor.n, source2.n)]

@named rc_model2 = System(rc_eqs2, t, systems = [resistor, capacitor, source2])
@test_throws ModelingToolkit.ExtraEquationsSystemException mtkcompile(rc_model2)
