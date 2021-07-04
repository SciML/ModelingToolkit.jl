using Test
using ModelingToolkit, OrdinaryDiffEq

@parameters t
@connector function Pin(;name)
    sts = @variables v(t)=1.0 i(t)=1.0
    ODESystem(Equation[], t, sts, []; name=name)
end

function ModelingToolkit.connect(::Type{Pin}, ps...)
    eqs = [
           0 ~ sum(p->p.i, ps) # KCL
          ]
    # KVL
    for i in 1:length(ps)-1
        push!(eqs, ps[i].v ~ ps[i+1].v)
    end

    return eqs
end

function Ground(;name)
    @named g = Pin()
    eqs = [g.v ~ 0]
    ODESystem(eqs, t, [], [], systems=[g]; name=name)
end

function OnePort(;name)
    @named p = Pin()
    @named n = Pin()
    sts = @variables v(t)=1.0 i(t)=1.0
    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           i ~ p.i
          ]
    ODESystem(eqs, t, sts, [], systems=[p, n]; name=name)
end

function ConstantVoltage(;name, V = 1.0)
    @named oneport = OnePort()
    @unpack v = oneport
    ps = @parameters V=V
    eqs = [
           V ~ v
          ]
    extend(oneport, ODESystem(eqs, t, [], ps; name=name); name=name)
end

function Resistor(;name, R = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters R=R
    eqs = [
           v ~ i * R
          ]
    extend(oneport, ODESystem(eqs, t, [], ps; name=name); name=name)
end

function Capacitor(;name, C = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters C=C
    D = Differential(t)
    eqs = [
           D(v) ~ i / C
          ]
    extend(oneport, ODESystem(eqs, t, [], ps; name=name); name=name)
end

function Inductor(; name, L = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters L=L
    D = Differential(t)
    eqs = [
           D(i) ~ v / L
          ]
    extend(oneport, ODESystem(eqs, t, [], ps; name=name); name=name)
end
