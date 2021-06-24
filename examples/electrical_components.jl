using Test
using ModelingToolkit, OrdinaryDiffEq

# Basic electric components
@parameters t
@connector function Pin(;name)
    @variables v(t)=1.0 i(t)=1.0
    ODESystem(Equation[], t, [v, i], [], name=name)
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
    ODESystem(eqs, t, [], [], systems=[g], name=name)
end

function Resistor(;name, R = 1.0)
    @named p = Pin()
    @named n = Pin()
    vars = @variables v(t)
    params = @parameters R=R
    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           v ~ p.i * R
          ]
    ODESystem(eqs, t, vars, params, systems=[p, n], name=name)
end

function Capacitor(; name, C = 1.0)
    @named p = Pin()
    @named n = Pin()
    vars = @variables v(t)
    params = @parameters C=C
    D = Differential(t)
    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           D(v) ~ p.i / C
          ]
    ODESystem(eqs, t, vars, params, systems=[p, n], name=name)
end

function ConstantVoltage(;name, V = 1.0)
    @named p = Pin()
    @named n = Pin()
    params = @parameters V=V
    eqs = [
           V ~ p.v - n.v
           0 ~ p.i + n.i
          ]
    ODESystem(eqs, t, [], params, systems=[p, n], name=name)
end

function Inductor(; name, L = 1.0)
    @named p = Pin()
    @named n = Pin()
    @variables v(t) i(t)
    @parameters L=L
    D = Differential(t)
    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           i ~ p.i
           D(i) ~ v / L
          ]
    ODESystem(eqs, t, [v, i], [L], systems=[p, n], name=name)
end
