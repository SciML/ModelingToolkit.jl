# Exposing More Parallelism By Tearing Algebraic Equations in ODESystems

Sometimes it can be very non-trivial to parallelize a system. In this tutorial,
we will demonstrate how to make use of `structural_simplify` to expose more
parallelism in the solution process and parallelize the resulting simulation.

## The Component Library

The following tutorial will use the following set of components describing
electrical circuits:

```@example tearing
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D

# Basic electric components
@connector function Pin(; name)
    @variables v(t)=1.0 i(t)=1.0 [connect = Flow]
    ODESystem(Equation[], t, [v, i], [], name = name)
end

function Ground(; name)
    @named g = Pin()
    eqs = [g.v ~ 0]
    compose(ODESystem(eqs, t, [], [], name = name), g)
end

function ConstantVoltage(; name, V = 1.0)
    val = V
    @named p = Pin()
    @named n = Pin()
    @parameters V = V
    eqs = [V ~ p.v - n.v
           0 ~ p.i + n.i]
    compose(ODESystem(eqs, t, [], [V], name = name), p, n)
end

@connector function HeatPort(; name)
    @variables T(t)=293.15 Q_flow(t)=0.0 [connect = Flow]
    ODESystem(Equation[], t, [T, Q_flow], [], name = name)
end

function HeatingResistor(; name, R = 1.0, TAmbient = 293.15, alpha = 1.0)
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

function HeatCapacitor(; name, rho = 8050, V = 1, cp = 460, TAmbient = 293.15)
    @parameters rho=rho V=V cp=cp
    C = rho * V * cp
    @named h = HeatPort()
    eqs = [
        D(h.T) ~ h.Q_flow / C
    ]
    compose(ODESystem(eqs, t, [], [rho, V, cp],
            name = name), h)
end

function Capacitor(; name, C = 1.0)
    @named p = Pin()
    @named n = Pin()
    @variables v(t) = 0.0
    @parameters C = C
    eqs = [v ~ p.v - n.v
           0 ~ p.i + n.i
           D(v) ~ p.i / C]
    compose(ODESystem(eqs, t, [v], [C],
            name = name), p, n)
end

function parallel_rc_model(i; name, source, ground, R, C)
    resistor = HeatingResistor(name = Symbol(:resistor, i), R = R)
    capacitor = Capacitor(name = Symbol(:capacitor, i), C = C)
    heat_capacitor = HeatCapacitor(name = Symbol(:heat_capacitor, i))

    rc_eqs = [connect(source.p, resistor.p)
              connect(resistor.n, capacitor.p)
              connect(capacitor.n, source.n, ground.g)
              connect(resistor.h, heat_capacitor.h)]

    compose(ODESystem(rc_eqs, t, name = Symbol(name, i)),
        [resistor, capacitor, source, ground, heat_capacitor])
end
```

## The Model

Assuming that the components are defined, our model is 50 resistors and
capacitors connected in parallel. Thus following the [acausal components tutorial](@ref acausal),
we can connect a bunch of RC components as follows:

```@example tearing
V = 2.0
@named source = ConstantVoltage(V = V)
@named ground = Ground()
N = 50
Rs = 10 .^ range(0, stop = -4, length = N)
Cs = 10 .^ range(-3, stop = 0, length = N)
rc_systems = map(1:N) do i
    parallel_rc_model(i; name = :rc, source = source, ground = ground, R = Rs[i], C = Cs[i])
end;
@variables E(t) = 0.0
eqs = [
    D(E) ~ sum(((i, sys),) -> getproperty(sys, Symbol(:resistor, i)).h.Q_flow,
    enumerate(rc_systems))
]
@named _big_rc = ODESystem(eqs, t, [E], [])
@named big_rc = compose(_big_rc, rc_systems)
```

Now let's say we want to expose a bit more parallelism via running tearing.
How do we do that?

```@example tearing
sys = structural_simplify(big_rc)
```

Done, that's it. There's no more to it.

## What Happened?

Yes, that's a good question! Let's investigate a little bit more what had happened.
If you look at the system we defined:

```@example tearing
length(equations(big_rc))
```

You see, it started as a massive 1051 set of equations. However, after eliminating
redundancies, we arrive at 151 equations:

```@example tearing
equations(sys)
```

That's not all though. In addition, the tearing process has turned the sets of
nonlinear equations into separate blocks and constructed a DAG for the dependencies
between the blocks. We can use the bipartite graph functionality to dig in and
investigate what this means:

```@example tearing
using ModelingToolkit.BipartiteGraphs
ts = TearingState(expand_connections(big_rc))
inc_org = BipartiteGraphs.incidence_matrix(ts.structure.graph)
blt_org = StructuralTransformations.sorted_incidence_matrix(ts, only_algeqs = true,
    only_algvars = true)
blt_reduced = StructuralTransformations.sorted_incidence_matrix(
    ModelingToolkit.get_tearing_state(sys),
    only_algeqs = true,
    only_algvars = true)
```

![](https://user-images.githubusercontent.com/1814174/110589027-d4ec9b00-8143-11eb-8880-651da986504d.PNG)

The figure on the left is the original incidence matrix of the algebraic equations.
Notice that the original formulation of the model has dependencies between different
equations, and so the full set of equations must be solved together. That exposes
no parallelism. However, the Block Lower Triangular (BLT) transformation exposes
independent blocks. This is then further improved by the tearing process, which
removes 90% of the equations and transforms the nonlinear equations into 50
independent blocks, *which can now all be solved in parallel*. The conclusion
is that, your attempts to parallelize are neigh: performing parallelism after
structural simplification greatly improves the problem that can be parallelized,
so this is better than trying to do it by hand.

After performing this, you can construct the `ODEProblem` and set
`parallel_form` to use the exposed parallelism in multithreaded function
constructions, but this showcases why `structural_simplify` is so important
to that process.
