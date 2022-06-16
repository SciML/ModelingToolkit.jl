using ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkit: PR_NULL_AFFECT, PeriodicEventCallback, PeriodicEventCallbacks, periodic_events


@variables t x(t)=1 v(t)=0 k(t)=0
D = Differential(t)


function affect!(u,p, t, state)
    p.R *= 2
    append!(state, p.R)
end

period = 1.0

@variables t
@connector function Pin(;name)
    sts = @variables v(t)=1.0 i(t)=1.0 [connect = Flow]
    ODESystem(Equation[], t, sts, []; name=name)
end

function Ground(;name)
    @named g = Pin()
    eqs = [g.v ~ 0]
    compose(ODESystem(eqs, t, [], []; name=name), g)
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
    compose(ODESystem(eqs, t, sts, []; name=name), p, n)
end

function Resistor(;name, R = 1.0, state=nothing)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters R=R
    eqs = [
           v ~ i * R
          ]
    extend(ODESystem(eqs, t, [], ps; name=name, periodic_events=(period, affect!, state)), oneport)
end

function Capacitor(;name, C = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters C=C
    D = Differential(t)
    eqs = [
           D(v) ~ i / C
          ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end

function ConstantVoltage(;name, V = 1.0)
    @named oneport = OnePort()
    @unpack v = oneport
    ps = @parameters V=V
    eqs = [
           V ~ v
          ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end

R = 1.0
C = 1.0
V = 1.0
mystate = Float64[]
@named resistor = Resistor(R=R; state=mystate)

@test length(ModelingToolkit.periodic_events(resistor)) == 1

@named capacitor = Capacitor(C=C)
@named source = ConstantVoltage(V=V)
@named ground = Ground()

rc_eqs = [
          connect(source.p, resistor.p)
          connect(resistor.n, capacitor.p)
          connect(capacitor.n, source.n)
          connect(capacitor.n, ground.g)
         ]

st2 = Float64[]

function affect2!(u,p, t, state)
    #@show(u, p, t, state)
end


@named _rc_model = ODESystem(rc_eqs, t, periodic_events=(5, affect2!))
@named rc_model = compose(_rc_model,
                          [resistor, capacitor, source, ground])

@test length(ModelingToolkit.periodic_events(rc_model)) == 2

sys = structural_simplify(rc_model)
@test length(ModelingToolkit.periodic_events(sys)) == 2
u0 = [
      capacitor.v => 0.0
     ]
prob = ODAEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Tsit5())
@test mystate == [2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0]
# plot(sol)

