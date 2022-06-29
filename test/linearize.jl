using ModelingToolkit

# r is an input, and y is an output.
@variables t x(t)=0 y(t)=0 u(t)=0 r(t)=0
@variables t x(t)=0 y(t)=0 u(t)=0 r(t)=0 [input = true]
@parameters kp = 1
D = Differential(t)

eqs = [u ~ kp * (r - y)
       D(x) ~ -x + u
       y ~ x]

@named sys = ODESystem(eqs, t)
linearize(sys, [r], [y])

##
```

  r ┌─────┐       ┌─────┐     ┌─────┐
───►│     ├──────►│     │  u  │     │
    │  F  │       │  C  ├────►│  P  │ y
    └─────┘     ┌►│     │     │     ├─┬─►
                │ └─────┘     └─────┘ │
                │                     │
                └─────────────────────┘
```

function plant(; name)
    @variables x(t) = 1
    @variables u(t)=0 [input = true] y(t)=0 [output = true]
    D = Differential(t)
    eqs = [D(x) ~ -x + u
           y ~ x]
    ODESystem(eqs, t; name = name)
end

function filt_(; name)
    @variables x(t)=0 y(t)=0 [output = true]
    @variables u(t)=0 [input = true]
    D = Differential(t)
    eqs = [D(x) ~ -2 * x + u
           y ~ x]
    ODESystem(eqs, t, name = name)
end

function controller(kp; name)
    @variables y(t)=0 r(t)=0 [input = true] u(t)=0
    @parameters kp = kp
    eqs = [
        u ~ kp * (r - y),
    ]
    ODESystem(eqs, t; name = name)
end

@named f = filt_()
@named c = controller(1)
@named p = plant()

connections = [f.y ~ c.r # filtered reference to controller reference
               c.u ~ p.u # controller output to plant input
               p.y ~ c.y]

@named cl = ODESystem(connections, t, systems = [f, c, p])

lin, xs = linearize(cl, cl.f.u, cl.p.x)

##
using ModelingToolkitStandardLibrary.Blocks: LimPID
#using ControlSystems
k = 400;
Ti = 0.5;
Td = 1;
Nd = 10;
#s = tf("s")
#expected_result_r = k*(1 + 1/(s*Ti)) |> ss
#expected_result_y = k*(1 + 1/(s*Ti) - s*Td / (1 + s*Td/N)) |> ss
@named pid = LimPID(; k, Ti, Td, Nd)
ModelingToolkit.unbound_inputs(pid)

@unpack reference, measurement, ctr_output = pid
lin = linearize(pid, [reference.u, measurement.u], [ctr_output.u])
lin, lin_fun = linearize(pid, [reference.u, measurement.u], [ctr_output.u]);
prob = ODEProblem(lin, [], (0.0, 1.0))
lin_fun(prob.u0, prob.p, 0.0)
