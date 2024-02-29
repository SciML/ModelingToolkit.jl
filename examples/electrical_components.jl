using Test
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D

@connector function Pin(; name)
    sts = @variables v(t) [guess = 1.0] i(t) [guess = 1.0, connect = Flow]
    ODESystem(Equation[], t, sts, []; name = name)
end

@component function Ground(; name)
    @named g = Pin()
    eqs = [g.v ~ 0]
    compose(ODESystem(eqs, t, [], []; name = name), g)
end

@component function OnePort(; name)
    @named p = Pin()
    @named n = Pin()
    sts = @variables v(t) [guess = 1.0] i(t) [guess = 1.0]
    eqs = [v ~ p.v - n.v
           0 ~ p.i + n.i
           i ~ p.i]
    compose(ODESystem(eqs, t, sts, []; name = name), p, n)
end

@component function Resistor(; name, R = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters R = R
    eqs = [
        v ~ i * R
    ]
    extend(ODESystem(eqs, t, [], ps; name = name), oneport)
end

@component function Capacitor(; name, C = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters C = C
    eqs = [
        D(v) ~ i / C
    ]
    extend(ODESystem(eqs, t, [], ps; name = name), oneport)
end

@component function ConstantVoltage(; name, V = 1.0)
    @named oneport = OnePort()
    @unpack v = oneport
    ps = @parameters V = V
    eqs = [
        V ~ v
    ]
    extend(ODESystem(eqs, t, [], ps; name = name), oneport)
end

@component function Inductor(; name, L = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters L = L
    eqs = [
        D(i) ~ v / L
    ]
    extend(ODESystem(eqs, t, [], ps; name = name), oneport)
end

@connector function HeatPort(; name)
    @variables T(t) [guess = 293.15] Q_flow(t) [guess = 0.0, connect = Flow]
    ODESystem(Equation[], t, [T, Q_flow], [], name = name)
end

@component function HeatingResistor(; name, R = 1.0, TAmbient = 293.15, alpha = 1.0)
    @named p = Pin()
    @named n = Pin()
    @named h = HeatPort()
    @variables v(t) RTherm(t)
    @parameters R=R TAmbient=TAmbient alpha=alpha
    eqs = [RTherm ~ R * (1 + alpha * (h.T - TAmbient))
           v ~ p.i * RTherm
           h.Q_flow ~ -v * p.i # -LossPower
           v ~ p.v - n.v
           0 ~ p.i + n.i]
    compose(ODESystem(eqs, t, [v, RTherm], [R, TAmbient, alpha],
            name = name), p, n, h)
end

@component function HeatCapacitor(; name, rho = 8050, V = 1, cp = 460, TAmbient = 293.15)
    @parameters rho=rho V=V cp=cp
    C = rho * V * cp
    @named h = HeatPort()
    eqs = [
        D(h.T) ~ h.Q_flow / C
    ]
    compose(ODESystem(eqs, t, [], [rho, V, cp],
            name = name), h)
end
