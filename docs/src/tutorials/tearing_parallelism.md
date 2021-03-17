# Exposing More Parallelism By Tearing Algebraic Equations in ODESystems

Sometimes it can be very non-trivial to parallelize a system. In this tutorial
we will demonstrate how to make use of `structural_simplify` to expose more
parallelism in the solution process and parallelize the resulting simulation.

## The Component Library

The following tutorial will use the following set of components describing
electrical circuits:

```julia
using ModelingToolkit, OrdinaryDiffEq

function connect_pin(ps...)
    eqs = [
           0 ~ sum(p->p.i, ps) # KCL
          ]
    # KVL
    for i in 1:length(ps)-1
        push!(eqs, ps[i].v ~ ps[i+1].v)
    end

    return eqs
end

function connect_heat(ps...)
    eqs = [
           0 ~ sum(p->p.Q_flow, ps) # KCL
          ]
    # KVL
    for i in 1:length(ps)-1
        push!(eqs, ps[i].T ~ ps[i+1].Q_flow)
    end

    return eqs
end

# Basic electric components
const t = Sym{ModelingToolkit.Parameter{Real}}(:t)
const D = Differential(t)
function Pin(;name)
    @variables v(t) i(t)
    ODESystem(Equation[], t, [v, i], [], name=name, defaults=Dict([v=>1.0, i=>1.0]))
end

function Ground(;name)
    @named g = Pin()
    eqs = [g.v ~ 0]
    ODESystem(eqs, t, [], [], systems=[g], name=name)
end

function ConstantVoltage(;name, V = 1.0)
    val = V
    @named p = Pin()
    @named n = Pin()
    @parameters V
    eqs = [
           V ~ p.v - n.v
           0 ~ p.i + n.i
          ]
    ODESystem(eqs, t, [], [V], systems=[p, n], defaults=Dict(V => val), name=name)
end

function HeatPort(;name)
    @variables T(t) Q_flow(t)
    return ODESystem(Equation[], t, [T, Q_flow], [], defaults=Dict(T=>293.15, Q_flow=>0.0), name=name)
end

function HeatingResistor(;name, R=1.0, TAmbient=293.15, alpha=1.0)
    R_val, TAmbient_val, alpha_val = R, TAmbient, alpha
    @named p = Pin()
    @named n = Pin()
    @named h = HeatPort()
    @variables v(t) RTherm(t)
    @parameters R TAmbient alpha
    eqs = [
           RTherm ~ R*(1 + alpha*(h.T - TAmbient))
           v ~ p.i * RTherm
           h.Q_flow ~ -v * p.i # -LossPower
           v ~ p.v - n.v
           0 ~ p.i + n.i
          ]
    ODESystem(
        eqs, t, [v, RTherm], [R, TAmbient, alpha], systems=[p, n, h],
        defaults=Dict(
            R=>R_val, TAmbient=>TAmbient_val, alpha=>alpha_val,
            v=>0.0, RTherm=>R_val
        ),
        name=name,
    )
end

function HeatCapacitor(;name, rho=8050, V=1, cp=460, TAmbient=293.15)
    rho_val, V_val, cp_val = rho, V, cp
    @parameters rho V cp
    C = rho*V*cp
    @named h = HeatPort()
    eqs = [
           D(h.T) ~ h.Q_flow / C
          ]
    ODESystem(
        eqs, t, [], [rho, V, cp], systems=[h],
        defaults=Dict(rho=>rho_val, V=>V_val, cp=>cp_val),
        name=name,
    )
end

function Capacitor(;name, C = 1.0)
    val = C
    @named p = Pin()
    @named n = Pin()
    @variables v(t)
    @parameters C
    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           D(v) ~ p.i / C
          ]
    ODESystem(
        eqs, t, [v], [C], systems=[p, n],
        defaults=Dict(v => 0.0, C => val),
        name=name
    )
end

function rc_model(i; name, source, ground, R, C)
    resistor = HeatingResistor(name=Symbol(:resistor, i), R=R)
    capacitor = Capacitor(name=Symbol(:capacitor, i), C=C)
    heat_capacitor = HeatCapacitor(name=Symbol(:heat_capacitor, i))

    rc_eqs = [
              connect_pin(source.p, resistor.p)
              connect_pin(resistor.n, capacitor.p)
              connect_pin(capacitor.n, source.n, ground.g)
              connect_heat(resistor.h, heat_capacitor.h)
             ]

    rc_model = ODESystem(rc_eqs, t, systems=[resistor, capacitor, source, ground, heat_capacitor], name=Symbol(name, i))
end
```

## The Model

Assuming that the components are defined, our model is 50 resistors and
capacitors connected in parallel. Thus following the [acausal components tutorial](@ref acausal),
we can connect a bunch of RC components as follows:

```julia
V = 2.0
source = ConstantVoltage(name=:source, V=V)
ground = Ground(name=:ground)
N = 50
Rs = 10 .^range(0, stop=-4, length=N)
Cs = 10 .^range(-3, stop=0, length=N)
rc_systems = map(1:N) do i
    rc_model(i; name=:rc, source=source, ground=ground, R=Rs[i], C=Cs[i])
end
@variables E(t)
eqs = [
       D(E) ~ sum(((i, sys),)->getproperty(sys, Symbol(:resistor, i)).h.Q_flow, enumerate(rc_systems))
      ]
big_rc = ODESystem(eqs, t, [], [], systems=rc_systems, defaults=Dict(E=>0.0))
```

Now let's say we want to expose a bit more parallelism via running tearing.
How do we do that?

```julia
sys = structural_simplify(big_rc)
```

Done, that's it. There's no more too it.

## What Happened?

Yes, that's a good question! Let's investigate a little bit more what had happened.
If you look at the system we defined:

```julia
equations(big_rc)

1051-element Vector{Equation}:
 Differential(t)(E(t)) ~ rc10₊resistor10₊h₊Q_flow(t) + rc11₊resistor11₊h₊Q_flow(t) + rc12₊resistor12₊h₊Q_flow(t) + rc13₊resistor13₊h₊Q_flow(t) + rc14₊resistor14₊h₊Q_flow(t) + rc15₊resistor15₊h₊Q_flow(t) + rc16₊resistor16₊h₊Q_flow(t) + rc17₊resistor17₊h₊Q_flow(t) + rc18₊resistor18₊h₊Q_flow(t) + rc19₊resistor19₊h₊Q_flow(t) + rc1₊resistor1₊h₊Q_flow(t) + rc20₊resistor20₊h₊Q_flow(t) + rc21₊resistor21₊h₊Q_flow(t) + rc22₊resistor22₊h₊Q_flow(t) + rc23₊resistor23₊h₊Q_flow(t) + rc24₊resistor24₊h₊Q_flow(t) + rc25₊resistor25₊h₊Q_flow(t) + rc26₊resistor26₊h₊Q_flow(t) + rc27₊resistor27₊h₊Q_flow(t) + rc28₊resistor28₊h₊Q_flow(t) + rc29₊resistor29₊h₊Q_flow(t) + rc2₊resistor2₊h₊Q_flow(t) + rc30₊resistor30₊h₊Q_flow(t) + rc31₊resistor31₊h₊Q_flow(t) + rc32₊resistor32₊h₊Q_flow(t) + rc33₊resistor33₊h₊Q_flow(t) + rc34₊resistor34₊h₊Q_flow(t) + rc35₊resistor35₊h₊Q_flow(t) + rc36₊resistor36₊h₊Q_flow(t) + rc37₊resistor37₊h₊Q_flow(t) + rc38₊resistor38₊h₊Q_flow(t) + rc39₊resistor39₊h₊Q_flow(t) + rc3₊resistor3₊h₊Q_flow(t) + rc40₊resistor40₊h₊Q_flow(t) + rc41₊resistor41₊h₊Q_flow(t) + rc42₊resistor42₊h₊Q_flow(t) + rc43₊resistor43₊h₊Q_flow(t) + rc44₊resistor44₊h₊Q_flow(t) + rc45₊resistor45₊h₊Q_flow(t) + rc46₊resistor46₊h₊Q_flow(t) + rc47₊resistor47₊h₊Q_flow(t) + rc48₊resistor48₊h₊Q_flow(t) + rc49₊resistor49₊h₊Q_flow(t) + rc4₊resistor4₊h₊Q_flow(t) + rc50₊resistor50₊h₊Q_flow(t) + rc5₊resistor5₊h₊Q_flow(t) + rc6₊resistor6₊h₊Q_flow(t) + rc7₊resistor7₊h₊Q_flow(t) + rc8₊resistor8₊h₊Q_flow(t) + rc9₊resistor9₊h₊Q_flow(t)
 0 ~ rc1₊resistor1₊p₊i(t) + rc1₊source₊p₊i(t)
 rc1₊source₊p₊v(t) ~ rc1₊resistor1₊p₊v(t)
 0 ~ rc1₊capacitor1₊p₊i(t) + rc1₊resistor1₊n₊i(t)
 ⋮
 rc50₊source₊V ~ rc50₊source₊p₊v(t) - (rc50₊source₊n₊v(t))
 0 ~ rc50₊source₊n₊i(t) + rc50₊source₊p₊i(t)
 rc50₊ground₊g₊v(t) ~ 0
 Differential(t)(rc50₊heat_capacitor50₊h₊T(t)) ~ rc50₊heat_capacitor50₊h₊Q_flow(t)*(rc50₊heat_capacitor50₊V^-1)*(rc50₊heat_capacitor50₊cp^-1)*(rc50₊heat_capacitor50₊rho^-1)
```

You see it started as a massive 1051 set of equations. However, after eliminating
redundancies we arrive at 151 equations:

```julia
equations(sys)

151-element Vector{Equation}:
 Differential(t)(E(t)) ~ rc1₊resistor1₊p₊i(t)*((rc1₊capacitor1₊v(t)) - rc1₊source₊V) + rc4₊resistor4₊p₊i(t)*((rc4₊capacitor4₊v(t)) - rc4₊source₊V) - ((rc10₊capacitor10₊p₊i(t))*(rc10₊source₊V - (rc10₊capacitor10₊v(t)))) - ((rc11₊capacitor11₊p₊i(t))*(rc11₊source₊V - (rc11₊capacitor11₊v(t)))) - ((rc12₊capacitor12₊p₊i(t))*(rc12₊source₊V - (rc12₊capacitor12₊v(t)))) - ((rc13₊capacitor13₊p₊i(t))*(rc13₊source₊V - (rc13₊capacitor13₊v(t)))) - ((rc14₊capacitor14₊p₊i(t))*(rc14₊source₊V - (rc14₊capacitor14₊v(t)))) - ((rc15₊capacitor15₊p₊i(t))*(rc15₊source₊V - (rc15₊capacitor15₊v(t)))) - ((rc16₊capacitor16₊p₊i(t))*(rc16₊source₊V - (rc16₊capacitor16₊v(t)))) - ((rc17₊capacitor17₊p₊i(t))*(rc17₊source₊V - (rc17₊capacitor17₊v(t)))) - ((rc18₊capacitor18₊p₊i(t))*(rc18₊source₊V - (rc18₊capacitor18₊v(t)))) - ((rc19₊capacitor19₊p₊i(t))*(rc19₊source₊V - (rc19₊capacitor19₊v(t)))) - ((rc20₊capacitor20₊p₊i(t))*(rc20₊source₊V - (rc20₊capacitor20₊v(t)))) - ((rc21₊capacitor21₊p₊i(t))*(rc21₊source₊V - (rc21₊capacitor21₊v(t)))) - ((rc22₊capacitor22₊p₊i(t))*(rc22₊source₊V - (rc22₊capacitor22₊v(t)))) - ((rc23₊capacitor23₊p₊i(t))*(rc23₊source₊V - (rc23₊capacitor23₊v(t)))) - ((rc24₊capacitor24₊p₊i(t))*(rc24₊source₊V - (rc24₊capacitor24₊v(t)))) - ((rc25₊capacitor25₊p₊i(t))*(rc25₊source₊V - (rc25₊capacitor25₊v(t)))) - ((rc26₊capacitor26₊p₊i(t))*(rc26₊source₊V - (rc26₊capacitor26₊v(t)))) - ((rc27₊capacitor27₊p₊i(t))*(rc27₊source₊V - (rc27₊capacitor27₊v(t)))) - ((rc28₊capacitor28₊p₊i(t))*(rc28₊source₊V - (rc28₊capacitor28₊v(t)))) - ((rc29₊capacitor29₊p₊i(t))*(rc29₊source₊V - (rc29₊capacitor29₊v(t)))) - ((rc2₊capacitor2₊p₊i(t))*(rc2₊source₊V - (rc2₊capacitor2₊v(t)))) - ((rc30₊capacitor30₊p₊i(t))*(rc30₊source₊V - (rc30₊capacitor30₊v(t)))) - ((rc31₊capacitor31₊p₊i(t))*(rc31₊source₊V - (rc31₊capacitor31₊v(t)))) - ((rc32₊capacitor32₊p₊i(t))*(rc32₊source₊V - (rc32₊capacitor32₊v(t)))) - ((rc33₊capacitor33₊p₊i(t))*(rc33₊source₊V - (rc33₊capacitor33₊v(t)))) - ((rc34₊capacitor34₊p₊i(t))*(rc34₊source₊V - (rc34₊capacitor34₊v(t)))) - ((rc35₊capacitor35₊p₊i(t))*(rc35₊source₊V - (rc35₊capacitor35₊v(t)))) - ((rc36₊capacitor36₊p₊i(t))*(rc36₊source₊V - (rc36₊capacitor36₊v(t)))) - ((rc37₊capacitor37₊p₊i(t))*(rc37₊source₊V - (rc37₊capacitor37₊v(t)))) - ((rc38₊capacitor38₊p₊i(t))*(rc38₊source₊V - (rc38₊capacitor38₊v(t)))) - ((rc39₊capacitor39₊p₊i(t))*(rc39₊source₊V - (rc39₊capacitor39₊v(t)))) - ((rc3₊capacitor3₊p₊i(t))*(rc3₊source₊V - (rc3₊capacitor3₊v(t)))) - ((rc40₊capacitor40₊p₊i(t))*(rc40₊source₊V - (rc40₊capacitor40₊v(t)))) - ((rc41₊capacitor41₊p₊i(t))*(rc41₊source₊V - (rc41₊capacitor41₊v(t)))) - ((rc42₊capacitor42₊p₊i(t))*(rc42₊source₊V - (rc42₊capacitor42₊v(t)))) - ((rc43₊capacitor43₊p₊i(t))*(rc43₊source₊V - (rc43₊capacitor43₊v(t)))) - ((rc44₊capacitor44₊p₊i(t))*(rc44₊source₊V - (rc44₊capacitor44₊v(t)))) - ((rc45₊capacitor45₊p₊i(t))*(rc45₊source₊V - (rc45₊capacitor45₊v(t)))) - ((rc46₊capacitor46₊p₊i(t))*(rc46₊source₊V - (rc46₊capacitor46₊v(t)))) - ((rc47₊capacitor47₊p₊i(t))*(rc47₊source₊V - (rc47₊capacitor47₊v(t)))) - ((rc48₊capacitor48₊p₊i(t))*(rc48₊source₊V - (rc48₊capacitor48₊v(t)))) - ((rc49₊capacitor49₊p₊i(t))*(rc49₊source₊V - (rc49₊capacitor49₊v(t)))) - ((rc50₊capacitor50₊p₊i(t))*(rc50₊source₊V - (rc50₊capacitor50₊v(t)))) - ((rc5₊capacitor5₊p₊i(t))*(rc5₊source₊V - (rc5₊capacitor5₊v(t)))) - ((rc6₊capacitor6₊p₊i(t))*(rc6₊source₊V - (rc6₊capacitor6₊v(t)))) - ((rc7₊capacitor7₊p₊i(t))*(rc7₊source₊V - (rc7₊capacitor7₊v(t)))) - ((rc8₊capacitor8₊p₊i(t))*(rc8₊source₊V - (rc8₊capacitor8₊v(t)))) - ((rc9₊capacitor9₊p₊i(t))*(rc9₊source₊V - (rc9₊capacitor9₊v(t))))
 0 ~ rc1₊resistor1₊R*rc1₊resistor1₊p₊i(t)*(1 + (rc1₊resistor1₊alpha*((-rc1₊resistor1₊TAmbient) - ((rc1₊resistor1₊p₊i(t))*((rc1₊capacitor1₊v(t)) - rc1₊source₊V))))) + rc1₊capacitor1₊v(t) - rc1₊source₊V
 Differential(t)(rc1₊capacitor1₊v(t)) ~ rc1₊resistor1₊p₊i(t)*(rc1₊capacitor1₊C^-1)
 Differential(t)(rc1₊heat_capacitor1₊h₊T(t)) ~ -rc1₊resistor1₊p₊i(t)*(rc1₊heat_capacitor1₊V^-1)*(rc1₊heat_capacitor1₊cp^-1)*(rc1₊heat_capacitor1₊rho^-1)*((rc1₊capacitor1₊v(t)) - rc1₊source₊V)
 ⋮
 Differential(t)(rc49₊heat_capacitor49₊h₊T(t)) ~ rc49₊capacitor49₊p₊i(t)*(rc49₊heat_capacitor49₊V^-1)*(rc49₊heat_capacitor49₊cp^-1)*(rc49₊heat_capacitor49₊rho^-1)*(rc49₊source₊V - (rc49₊capacitor49₊v(t)))
 0 ~ rc50₊resistor50₊R*rc50₊capacitor50₊p₊i(t)*(1 + (rc50₊resistor50₊alpha*(((rc50₊capacitor50₊p₊i(t))*(rc50₊source₊V - (rc50₊capacitor50₊v(t)))) - rc50₊resistor50₊TAmbient))) - (rc50₊source₊V - (rc50₊capacitor50₊v(t)))
 Differential(t)(rc50₊capacitor50₊v(t)) ~ rc50₊capacitor50₊p₊i(t)*(rc50₊capacitor50₊C^-1)
 Differential(t)(rc50₊heat_capacitor50₊h₊T(t)) ~ rc50₊capacitor50₊p₊i(t)*(rc50₊heat_capacitor50₊V^-1)*(rc50₊heat_capacitor50₊cp^-1)*(rc50₊heat_capacitor50₊rho^-1)*(rc50₊source₊V - (rc50₊capacitor50₊v(t)))
```

That's not all though. In addition, the tearing process has turned the sets of
nonlinear equations into separate blocks and constructed a DAG for the dependencies
between the blocks. We can use the bipartite graph functionality to dig in and
investigate what this means:

```julia
using ModelingToolkit.BipartiteGraphs
big_rc = initialize_system_structure(big_rc)
inc_org = BipartiteGraphs.incidence_matrix(structure(big_rc).graph)
blt_org = StructuralTransformations.sorted_incidence_matrix(big_rc, only_algeqs=true, only_algvars=true)
blt_reduced = StructuralTransformations.sorted_incidence_matrix(sys, only_algeqs=true, only_algvars=true)
```

![](https://user-images.githubusercontent.com/1814174/110589027-d4ec9b00-8143-11eb-8880-651da986504d.PNG)

The figure on the left is the original incidence matrix of the algebraic equations.
Notice that the original formulation of the model has dependencies between different
equations, and so the full set of equations must be solved together. That exposes
no parallelism. However, the Block Lower Triangular (BLT) transformation exposes
independent blocks. This is then further impoved by the tearing process, which
removes 90% of the equations and transforms the nonlinear equations into 50
independent blocks *which can now all be solved in parallel*. The conclusion
is that, your attempts to parallelize are neigh: performing parallelism after
structural simplification greatly improves the problem that can be parallelized,
so this is better than trying to do it by hand.

After performing this, you can construct the `ODEProblem`/`ODAEProblem` and set
`parallel_form` to use the exposed parallelism in multithreaded function
constructions, but this showcases why `structural_simplify` is so important
to that process.
