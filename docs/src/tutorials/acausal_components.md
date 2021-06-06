# Acausal Component-Based Modeling the RC Circuit

In this tutorial we will build a hierarchical acausal component-based model of
the RC circuit. The RC circuit is a simple example where we connect a resistor
and a capacitor. [Kirchoff's laws](https://en.wikipedia.org/wiki/Kirchhoff%27s_circuit_laws)
are then applied to state equalities between currents and voltages. This
specifies a differential-algebraic equation (DAE) system, where the algebraic
equations are given by the constraints and equalities between different
component variables. We then simplify this to an ODE by eliminating the
equalities before solving. Let's see this in action.

## Copy-Paste Example

```julia
using ModelingToolkit, Plots, DifferentialEquations

@parameters t

# Basic electric components
function Pin(;name)
    @variables v(t) i(t)
    ODESystem(Equation[], t, [v, i], [], name=name, defaults=[v=>1.0, i=>1.0])
end

function Ground(;name)
    @named g = Pin()
    eqs = [g.v ~ 0]
    ODESystem(eqs, t, [], [], systems=[g], name=name)
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

function Capacitor(; name, C = 1.0)
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

R = 1.0
C = 1.0
V = 1.0
@named resistor = Resistor(R=R)
@named capacitor = Capacitor(C=C)
@named source = ConstantVoltage(V=V)
@named ground = Ground()

function connect_pins(ps...)
    eqs = [
           0 ~ sum(p->p.i, ps) # KCL
          ]
    # KVL
    for i in 1:length(ps)-1
        push!(eqs, ps[i].v ~ ps[i+1].v)
    end

    return eqs
end

rc_eqs = [
          connect_pins(source.p, resistor.p)
          connect_pins(resistor.n, capacitor.p)
          connect_pins(capacitor.n, source.n, ground.g)
         ]

@named rc_model = ODESystem(rc_eqs, t,
                            systems=[resistor, capacitor, source, ground])
sys = structural_simplify(rc_model)
u0 = [
      capacitor.v => 0.0
      capacitor.p.i => 0.0
     ]
prob = ODAEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Tsit5())
plot(sol)
```

![](https://user-images.githubusercontent.com/1814174/109416294-55184100-798b-11eb-9f05-766a793f0ba2.png)

## Explanation

### Building the Component Library

For each of our components we use a Julia function which emits an `ODESystem`.
At the top we start with defining the fundamental qualities of an electrical
circuit component. At every input and output pin a circuit component has
two values: the current at the pin and the voltage. Thus we define the `Pin`
component to simply be the values there:

```julia
function Pin(;name)
    @variables v(t) i(t)
    ODESystem(Equation[], t, [v, i], [], name=name, defaults=[v=>1.0, i=>1.0])
end
```

Note that this is an incompletely specified ODESystem: it cannot be simulated
on its own because the equations for `v(t)` and `i(t)` are unknown. Instead
this just gives a common syntax for receiving this pair with some default
values. Notice that in a component we define the `name` as a keyword argument:
this is because later we will generate different `Pin` objects with different
names to correspond to duplicates of this topology with unique variables.
One can then construct a `Pin` like:

```julia
Pin(name=:mypin1)
```

or equivalently using the `@named` helper macro:

```julia
@named mypin1 = Pin()
```

Next we build our ground node. A ground node is just a pin that is connected
to a constant voltage reservoir, typically taken to be `V=0`. Thus to define
this component, we generate an `ODESystem` with a `Pin` subcomponent and specify
that the voltage in such a `Pin` is equal to zero. This gives:

```julia
function Ground(;name)
    @named g = Pin()
    eqs = [g.v ~ 0]
    ODESystem(eqs, t, [], [], systems=[g], name=name)
end
```

Next we build a resistor. A resistor is an object that has two `Pin`s, the positive
and the negative pins, and follows Ohm's law: `v = i*r`. The voltage of the
resistor is given as the voltage difference across the two pins while by conservation
of charge we know that the current in must equal the current out, which means
(no matter the direction of the current flow) the sum of the currents must be
zero. This leads to our resistor equations:

```julia
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
```

Notice that we have created this system with a `defaults` for the resistor's
resistance. By doing so, if the resistance of this resistor is not overridden
by a higher level default or overridden at `ODEProblem` construction time, this
will be the value of the resistance.

Using our knowledge of circuits we similarly construct the Capacitor:

```julia
function Capacitor(; name, C = 1.0)
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
```

Now we want to build a constant voltage electrical source term. We can think of
this as similarly being a two pin object, where the object itself is kept at a
constant voltage, essentially generating the electrical current. We would then
model this as:

```julia
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
```

### Connecting and Simulating Our Electrical Circuit

Now we are ready to simulate our circuit. Let's build our four components:
a `resistor`, `capacitor`, `source`, and `ground` term. For simplicity we will
make all of our parameter values 1. This is done by:

```julia
R = 1.0
C = 1.0
V = 1.0
@named resistor = Resistor(R=R)
@named capacitor = Capacitor(C=C)
@named source = ConstantVoltage(V=V)
@named ground = Ground()
```

Next we have to define how we connect the circuit. Whenever two `Pin`s in a
circuit are connected together, the system satisfies
[Kirchoff's laws](https://en.wikipedia.org/wiki/Kirchhoff%27s_circuit_laws),
i.e. that currents sum to zero and voltages across the pins are equal. Thus
we will build a helper function `connect_pins` which implements these rules:

```julia
function connect_pins(ps...)
    eqs = [
           0 ~ sum(p->p.i, ps) # KCL
          ]
    # KVL
    for i in 1:length(ps)-1
        push!(eqs, ps[i].v ~ ps[i+1].v)
    end

    return eqs
end
```

Finally we will connect the pieces of our circuit together. Let's connect the
positive pin of the resistor to the source, the negative pin of the resistor
to the capacitor, and the negative pin of the capacitor to a junction between
the source and the ground. This would mean our connection equations are:

```julia
rc_eqs = [
          connect_pins(source.p, resistor.p)
          connect_pins(resistor.n, capacitor.p)
          connect_pins(capacitor.n, source.n, ground.g)
         ]
```

Finally we build our four component model with these connection rules:

```julia
@named rc_model = ODESystem(rc_eqs, t,
                            systems=[resistor, capacitor, source, ground])
```

Notice that this model is acasual because we have not specified anything about
the causality of the model. We have simply specified what is true about each
of the variables. This forms a system of differential-algebraic equations
(DAEs) which define the evolution of each state of the system. The
equations are:

```julia
equations(rc_model)

16-element Vector{Equation}:
 0 ~ resistor₊p₊i(t) + source₊p₊i(t)
 source₊p₊v(t) ~ resistor₊p₊v(t)
 0 ~ capacitor₊p₊i(t) + resistor₊n₊i(t)
 resistor₊n₊v(t) ~ capacitor₊p₊v(t)
 ⋮
 Differential(t)(capacitor₊v(t)) ~ capacitor₊p₊i(t)*(capacitor₊C^-1)
 source₊V ~ source₊p₊v(t) - (source₊n₊v(t))
 0 ~ source₊n₊i(t) + source₊p₊i(t)
 ground₊g₊v(t) ~ 0
```

the states are:

```julia
states(rc_model)

16-element Vector{Term{Real}}:
 resistor₊p₊i(t)
 source₊p₊i(t)
 source₊p₊v(t)
 resistor₊p₊v(t)
 ⋮
 source₊n₊v(t)
 ground₊g₊v(t)
 resistor₊v(t)
 capacitor₊v(t)
```

and the parameters are:

```julia
parameters(rc_model)

3-element Vector{Any}:
 resistor₊R
 capacitor₊C
 source₊V
```

## Simplifying and Solving this System

This system could be solved directly as a DAE using [one of the DAE solvers
from DifferentialEquations.jl](https://diffeq.sciml.ai/stable/solvers/dae_solve/).
However, let's take a second to symbolically simplify the system before doing the
solve. The function `structural_simplify` looks for all of the equalities and
eliminates unnecessary variables to build the leanest numerical representation
of the system. Let's see what it does here:

```julia
sys = structural_simplify(rc_model)
equations(sys)

2-element Vector{Equation}:
 0 ~ capacitor₊v(t) + resistor₊R*capacitor₊p₊i(t) - source₊V
 Differential(t)(capacitor₊v(t)) ~ capacitor₊p₊i(t)*(capacitor₊C^-1)
```

```julia
states(sys)

2-element Vector{Any}:
 capacitor₊v(t)
 capacitor₊p₊i(t)
```

After structural simplification we are left with a system of only two equations
with two state variables. One of the equations is a differential equation
while the other is an algebraic equation. We can then give the values for the
initial conditions of our states and solve the system by converting it to
an ODEProblem in mass matrix form and solving it with an [ODEProblem mass matrix
DAE solver](https://diffeq.sciml.ai/stable/solvers/dae_solve/#OrdinaryDiffEq.jl-(Mass-Matrix)).
This is done as follows:

```julia
u0 = [
      capacitor.v => 0.0
      capacitor.p.i => 0.0
     ]
prob = ODEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Rodas4())
plot(sol)
```

![](https://user-images.githubusercontent.com/1814174/109416295-55184100-798b-11eb-96d1-5bb7e40135ba.png)

However, we can also choose to use the "torn nonlinear system" to remove all
of the algebraic variables from the solution of the system. Note that this
requires having done `structural_simplify`. This is done by using `ODAEProblem`
like:

```julia
u0 = [
      capacitor.v => 0.0
      capacitor.p.i => 0.0
     ]
prob = ODAEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Rodas4())
plot(sol)
```

![](https://user-images.githubusercontent.com/1814174/109416294-55184100-798b-11eb-9f05-766a793f0ba2.png)

Notice that this solves the whole system by only solving for one variable!

However, what if we wanted to plot the timeseries of a different variable? Do
not worry, that information was not thrown away! Instead, transformations
like `structural_simplify` simply change state variables into `observed`
variables. Let's see what our observed variables are:

```julia
observed(sys)

14-element Vector{Equation}:
 resistor₊p₊i(t) ~ capacitor₊p₊i(t)
 capacitor₊n₊v(t) ~ 0.0
 source₊n₊v(t) ~ 0.0
 ground₊g₊i(t) ~ 0.0
 source₊n₊i(t) ~ capacitor₊p₊i(t)
 source₊p₊i(t) ~ -capacitor₊p₊i(t)
 capacitor₊n₊i(t) ~ -capacitor₊p₊i(t)
 resistor₊n₊i(t) ~ -capacitor₊p₊i(t)
 ground₊g₊v(t) ~ 0.0
 source₊p₊v(t) ~ source₊V
 capacitor₊p₊v(t) ~ capacitor₊v(t)
 resistor₊p₊v(t) ~ source₊p₊v(t)
 resistor₊n₊v(t) ~ capacitor₊p₊v(t)
 resistor₊v(t) ~ -((capacitor₊p₊v(t)) - (source₊p₊v(t)))
```

These are explicit algebraic equations which can then be used to reconstruct
the required variables on the fly. This leads to dramatic computational savings
because implicitly solving an ODE scales like O(n^3), so making there be as
few states as possible is good!

The solution object can be accessed via its symbols. For example, let's retrieve
the voltage of the resistor over time:

```julia
sol[resistor.v]
```

or we can plot the timeseries of the resistor's voltage:

```julia
plot(sol, vars=[resistor.v])
```
