using Test
using ModelingToolkit#, OrdinaryDiffEq
import ModelingToolkit: ExtraVariablesSystemException, ExtraEquationsSystemException

@parameters t
function Pin(;name)
    @variables v(t) i(t)
    ODESystem(Equation[], t, [v, i], [], name=name, defaults=[v=>1.0, i=>1.0])
end

function UnderdefinedConstantVoltage(;name, V = 1.0)
    val = V
    @named p = Pin()
    @named n = Pin()
    @parameters V
    eqs = [
           V ~ p.v - n.v
           #0 ~ p.i + n.i
          ]
    ODESystem(eqs, t, [], [V], systems=[p, n], defaults=Dict(V => val), name=name)
end

function OverdefinedConstantVoltage(;name, V = 1.0, I = 1.0)
    val = V
    val2 = I
    @named p = Pin()
    @named n = Pin()
    @parameters V I
    eqs = [
           V ~ p.v - n.v
           n.i ~ I
           p.i ~ I
          ]
    ODESystem(eqs, t, [], [V], systems=[p, n], defaults=Dict(V => val, I => val2), name=name)
end

function Resistor(;name, R = 1.0)
    val = R
    @named p = Pin()
    @named n = Pin()
    @variables v(t)
    @parameters R
    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           v ~ p.i * R
          ]
    ODESystem(eqs, t, [v], [R], systems=[p, n], defaults=Dict(R => val), name=name)
end

function Capacitor(;name, C = 1.0)
    val = C
    @named p = Pin()
    @named n = Pin()
    @variables v(t)
    @parameters C
    D = Differential(t)
    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           D(v) ~ p.i / C
          ]
    ODESystem(eqs, t, [v], [C], systems=[p, n], defaults=Dict(C => val), name=name)
end

function ModelingToolkit.connect(ps...)
    eqs = [
           0 ~ sum(p->p.i, ps) # KCL
          ]
    # KVL
    for i in 1:length(ps)-1
        push!(eqs, ps[i].v ~ ps[i+1].v)
    end

    return eqs
end

R = 1.0
C = 1.0
V = 1.0
@named resistor = Resistor(R=R)
@named capacitor = Capacitor(C=C)
@named source = UnderdefinedConstantVoltage(V=V)

rc_eqs = [
          connect(source.p, resistor.p)
          connect(resistor.n, capacitor.p)
          connect(capacitor.n, source.n)
         ]

@named rc_model = ODESystem(rc_eqs, t, systems=[resistor, capacitor, source])


@test_throws ModelingToolkit.ExtraVariablesSystemException structural_simplify(rc_model)


@named source2 = OverdefinedConstantVoltage(V=V, I=V/R)
rc_eqs2 = [
          connect(source2.p, resistor.p)
          connect(resistor.n, capacitor.p)
          connect(capacitor.n, source2.n)
         ]

@named rc_model2 = ODESystem(rc_eqs2, t, systems=[resistor, capacitor, source2])
@test_throws ModelingToolkit.ExtraEquationsSystemException structural_simplify(rc_model2)
