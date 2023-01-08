# [Acausal Component-Based Modeling](@id acausal)

In this tutorial, we will build a hierarchical acausal component-based model of
the RC circuit. The RC circuit is a simple example where we connect a resistor
and a capacitor. [Kirchhoff's laws](https://en.wikipedia.org/wiki/Kirchhoff%27s_circuit_laws)
are then applied to state equalities between currents and voltages. This
specifies a differential-algebraic equation (DAE) system, where the algebraic
equations are given by the constraints and equalities between different
component variables. We then simplify this to an ODE by eliminating the
equalities before solving. Let's see this in action.

!!! note

    This tutorial teaches how to build the entire RC circuit from scratch.
    However, to simulate electric components with more ease, check out the
    [ModelingToolkitStandardLibrary.jl](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/)
    which includes a
    [tutorial for simulating RC circuits with pre-built components](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/tutorials/rc_circuit/)

## Copy-Paste Example

```@example acausal
using ModelingToolkit, Plots, DifferentialEquations

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

function Resistor(;name, R = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters R=R
    eqs = [
           v ~ i * R
          ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
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
@named resistor = Resistor(R=R)
@named capacitor = Capacitor(C=C)
@named source = ConstantVoltage(V=V)
@named ground = Ground()

rc_eqs = [
          connect(source.p, resistor.p)
          connect(resistor.n, capacitor.p)
          connect(capacitor.n, source.n)
          connect(capacitor.n, ground.g)
         ]

@named _rc_model = ODESystem(rc_eqs, t)
@named rc_model = compose(_rc_model,
                          [resistor, capacitor, source, ground])
sys = structural_simplify(rc_model)
u0 = [
      capacitor.v => 0.0
     ]
prob = ODAEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Tsit5())
plot(sol)
```

## Explanation

We wish to build the following RC circuit by building individual components and connecting the pins:

![](https://user-images.githubusercontent.com/1814174/172466302-907d39f3-6d2c-4d16-84a8-6de32bca757e.png)

### Building the Component Library

For each of our components, we use a Julia function which emits an `ODESystem`.
At the top, we start with defining the fundamental qualities of an electric
circuit component. At every input and output pin, a circuit component has
two values: the current at the pin and the voltage. Thus we define the `Pin`
component (connector) to simply be the values there. Whenever two `Pin`s in a
circuit are connected together, the system satisfies [Kirchhoff's laws](https: //en.wikipedia.org/wiki/Kirchhoff%27s_circuit_laws),
i.e. that currents sum to zero and voltages across the pins are equal.
`[connect = Flow]` informs MTK that currents ought to sum to zero, and by
default, variables are equal in a connection.

```@example acausal
@connector function Pin(;name)
    sts = @variables v(t)=1.0 i(t)=1.0 [connect = Flow]
    ODESystem(Equation[], t, sts, []; name=name)
end
```

Note that this is an incompletely specified ODESystem: it cannot be simulated
on its own because the equations for `v(t)` and `i(t)` are unknown. Instead,
this just gives a common syntax for receiving this pair with some default
values. Notice that in a component, we define the `name` as a keyword argument:
this is because later we will generate different `Pin` objects with different
names to correspond to duplicates of this topology with unique variables.
One can then construct a `Pin` like:

```@example acausal
Pin(name=:mypin1)
```

or equivalently, using the `@named` helper macro:

```@example acausal
@named mypin1 = Pin()
```

Next, we build our ground node. A ground node is just a pin that is connected
to a constant voltage reservoir, typically taken to be `V=0`. Thus to define
this component, we generate an `ODESystem` with a `Pin` subcomponent and specify
that the voltage in such a `Pin` is equal to zero. This gives:

```@example acausal
function Ground(;name)
    @named g = Pin()
    eqs = [g.v ~ 0]
    compose(ODESystem(eqs, t, [], []; name=name), g)
end
```

Next we build a `OnePort`: an abstraction for all simple electric component
with two pins. The voltage difference between the positive pin and the negative
pin is the voltage of the component, the current between two pins must sum to
zero, and the current of the component equals to the current of the positive
pin.

```@example acausal
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
```

Next we build a resistor. A resistor is an object that has two `Pin`s, the positive
and the negative pins, and follows Ohm's law: `v = i*r`. The voltage of the
resistor is given as the voltage difference across the two pins, while by conservation
of charge we know that the current in must equal the current out, which means
(no matter the direction of the current flow) the sum of the currents must be
zero. This leads to our resistor equations:

```@example acausal
function Resistor(;name, R = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters R=R
    eqs = [
           v ~ i * R
          ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end
```

Notice that we have created this system with a default parameter `R` for the
resistor's resistance. By doing so, if the resistance of this resistor is not
overridden by a higher level default or overridden at `ODEProblem` construction
time, this will be the value of the resistance. Also, note the use of `@unpack`
and `extend`. For the `Resistor`, we want to simply inherit `OnePort`'s
equations and states and extend them with a new equation. ModelingToolkit makes
a new namespaced variable `oneport₊v(t)` when using the syntax `oneport.v`, and
we can use `@unpack` to avoid the namespacing.

Using our knowledge of circuits, we similarly construct the `Capacitor`:

```@example acausal
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
```

Now we want to build a constant voltage electric source term. We can think of
this as similarly being a two pin object, where the object itself is kept at a
constant voltage, essentially generating the electric current. We would then
model this as:

```@example acausal
function ConstantVoltage(;name, V = 1.0)
    @named oneport = OnePort()
    @unpack v = oneport
    ps = @parameters V=V
    eqs = [
           V ~ v
          ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end
```

### Connecting and Simulating Our Electric Circuit

Now we are ready to simulate our circuit. Let's build our four components:
a `resistor`, `capacitor`, `source`, and `ground` term. For simplicity, we will
make all of our parameter values 1. This is done by:

```@example acausal
R = 1.0
C = 1.0
V = 1.0
@named resistor = Resistor(R=R)
@named capacitor = Capacitor(C=C)
@named source = ConstantVoltage(V=V)
@named ground = Ground()
```

Subsequently, we will connect the pieces of our circuit together. Let's connect the
positive pin of the resistor to the source, the negative pin of the resistor
to the capacitor, and the negative pin of the capacitor to a junction between
the source and the ground. This would mean our connection equations are:

```@example acausal
rc_eqs = [
          connect(source.p, resistor.p)
          connect(resistor.n, capacitor.p)
          connect(capacitor.n, source.n)
          connect(capacitor.n, ground.g)
         ]
```

Finally, we build our four-component model with these connection rules:

```@example acausal
@named _rc_model = ODESystem(rc_eqs, t)
@named rc_model = compose(_rc_model,
                          [resistor, capacitor, source, ground])
```

Note that we can also specify the subsystems in a vector. This model is acausal
because we have not specified anything about the causality of the model. We have
simply specified what is true about each of the variables. This forms a system
of differential-algebraic equations (DAEs) which define the evolution of each
state of the system. The equations are:

```@example acausal
equations(expand_connections(rc_model))
```

the states are:

```@example acausal
states(rc_model)
```

and the parameters are:

```@example acausal
parameters(rc_model)
```

## Simplifying and Solving this System

This system could be solved directly as a DAE using [one of the DAE solvers
from DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/solvers/dae_solve/).
However, let's take a second to symbolically simplify the system before doing the
solve. Although we can use ODE solvers that handle mass matrices to solve the
above system directly, we want to run the `structural_simplify` function first,
as it eliminates many unnecessary variables to build the leanest numerical
representation of the system. Let's see what it does here:

```@example acausal
sys = structural_simplify(rc_model)
equations(sys)
```

```@example acausal
states(sys)
```

After structural simplification, we are left with a system of only two equations
with two state variables. One of the equations is a differential equation,
while the other is an algebraic equation. We can then give the values for the
initial conditions of our states, and solve the system by converting it to
an ODEProblem in mass matrix form and solving it with an [ODEProblem mass matrix
DAE solver](https://docs.sciml.ai/DiffEqDocs/stable/solvers/dae_solve/#OrdinaryDiffEq.jl-(Mass-Matrix)).
This is done as follows:

```@example acausal
u0 = [
      capacitor.v => 0.0
      capacitor.p.i => 0.0
     ]
prob = ODEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Rodas4())
plot(sol)
```

Since we have run `structural_simplify`, MTK can numerically solve all the
unreduced algebraic equations using the `ODAEProblem` (note the
letter `A`):

```@example acausal
u0 = [
      capacitor.v => 0.0
     ]
prob = ODAEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Rodas4())
plot(sol)
```

Notice that this solves the whole system by only solving for one variable!

However, what if we wanted to plot the timeseries of a different variable? Do
not worry, that information was not thrown away! Instead, transformations
like `structural_simplify` simply change state variables into `observed`
variables. Let's see what our observed variables are:

```@example acausal
observed(sys)
```

These are explicit algebraic equations which can then be used to reconstruct
the required variables on the fly. This leads to dramatic computational savings
because implicitly solving an ODE scales like O(n^3), so making there be as
few states as possible is good!

The solution object can be accessed via its symbols. For example, let's retrieve
the voltage of the resistor over time:

```@example acausal
sol[resistor.v]
```

or we can plot the timeseries of the resistor's voltage:

```@example acausal
plot(sol, idxs=[resistor.v])
```
