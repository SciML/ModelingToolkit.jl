using Test
using ModelingToolkit, OrdinaryDiffEq

# Basic electric components
@parameters t
@connector function Pin(;name)
    sts = @variables v(t)=1.0 i(t)=1.0
    ODESystem(Equation[], t, sts, [], name=name)
end

@namespace function ModelingToolkit.connect(::Type{Pin}, ps...)
    eqs = [
           0 ~ sum(p->p.i, ps) # KCL
          ]
    # KVL
    for i in 1:length(ps)-1
        push!(eqs, ps[i].v ~ ps[i+1].v)
    end

    return eqs
end

@namespace function Ground(;name)
    @named g = Pin()
    eqs = [g.v ~ 0]
    ODESystem(eqs, t, [], [], systems=[g], name=name)
end

@namespace function ConstantVoltage(;name, V = 1.0)
    @named p = Pin()
    @named n = Pin()
    ps = @parameters V=V
    eqs = [
           V ~ p.v - n.v
           0 ~ p.i + n.i
          ]
    ODESystem(eqs, t, [], ps, systems=[p, n], name=name)
end

@namespace function Resistor(;name, R = 1.0)
    @named p = Pin()
    @named n = Pin()
    sts = @variables v(t)
    ps = @parameters R=R
    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           v ~ p.i * R
          ]
    ODESystem(eqs, t, sts, ps, systems=[p, n], name=name)
end

@namespace function Capacitor(;name, C = 1.0)
    @named p = Pin()
    @named n = Pin()
    sts = @variables v(t)
    ps = @parameters C=C
    D = Differential(t)
    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           D(v) ~ p.i / C
          ]
    ODESystem(eqs, t, sts, ps, systems=[p, n], name=name)
end

@namespace function Inductor(; name, L = 1.0)
    @named p = Pin()
    @named n = Pin()
    sts = @variables v(t) i(t)
    ps = @parameters L=L
    D = Differential(t)
    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           i ~ p.i
           D(i) ~ v / L
          ]
    ODESystem(eqs, t, sts, ps, systems=[p, n], name=name)
end
