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
using ModelingToolkit, Plots, OrdinaryDiffEq, Setfield
using ModelingToolkit: t_nounits as t, D_nounits as D

# Define the Pin connector
function Pin(; name)
    @variables v(t) i(t) [connect = Flow]
    sys = System(Equation[], t, [v, i], []; name)
    sys = @set sys.connector_type = ModelingToolkit.connector_type(sys)
    return sys
end

# Define Ground component
function Ground(; name)
    @named g = Pin()
    eqs = [g.v ~ 0]
    System(eqs, t, [], []; systems = [g], name)
end

# Define OnePort base component
function OnePort(; name)
    @named p = Pin()
    @named n = Pin()
    @variables v(t) i(t)
    eqs = [
        v ~ p.v - n.v
        0 ~ p.i + n.i
        i ~ p.i
    ]
    System(eqs, t, [v, i], []; systems = [p, n], name)
end

# Define Resistor component
function Resistor(; name, R = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    @parameters R = R # Sets the default resistance
    eqs = [v ~ i * R]
    extend(System(eqs, t, [], [R]; name), oneport)
end

# Define Capacitor component
function Capacitor(; name, C = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    @parameters C = C
    eqs = [D(v) ~ i / C]
    extend(System(eqs, t, [], [C]; name), oneport)
end

# Define ConstantVoltage source
function ConstantVoltage(; name, V = 1.0)
    @named oneport = OnePort()
    @unpack v = oneport
    @parameters V = V
    eqs = [V ~ v]
    extend(System(eqs, t, [], [V]; name), oneport)
end

# Build the RC circuit
@named resistor = Resistor(R = 2.0)
@named capacitor = Capacitor(C = 1.0)
@named source = ConstantVoltage(V = 1.0)
@named ground = Ground()

rc_eqs = [
    connect(source.p, resistor.p)
    connect(resistor.n, capacitor.p)
    connect(capacitor.n, source.n)
    connect(capacitor.n, ground.g)
]

@named rc_model = System(rc_eqs, t; systems = [resistor, capacitor, source, ground])
rc_model = mtkcompile(rc_model)

u0 = [
    rc_model.capacitor.v => 0.0
]
prob = ODEProblem(rc_model, u0, (0, 10.0))
sol = solve(prob)
plot(sol)
```

## Explanation

We wish to build the following RC circuit by building individual components and connecting the pins:

![](https://user-images.githubusercontent.com/1814174/172466302-907d39f3-6d2c-4d16-84a8-6de32bca757e.png)

### Building the Component Library

For each of our components, we define a function that returns a `System`.
At the top, we start with defining the fundamental qualities of an electric
circuit component. At every input and output pin, a circuit component has
two values: the current at the pin and the voltage. Thus we define the `Pin`
connector to simply be the values there. Whenever two `Pin`s in a
circuit are connected together, the system satisfies [Kirchhoff's laws](https://en.wikipedia.org/wiki/Kirchhoff%27s_circuit_laws),
i.e. that currents sum to zero and voltages across the pins are equal.
`[connect = Flow]` informs MTK that currents ought to sum to zero, and by
default, variables are equal in a connection.

```@example acausal
function Pin(; name)
    @variables v(t) i(t) [connect = Flow]
    sys = System(Equation[], t, [v, i], []; name)
    sys = @set sys.connector_type = ModelingToolkit.connector_type(sys)
    return sys
end
```

Note that this is an incompletely specified System: it cannot be simulated
on its own because the equations for `v(t)` and `i(t)` are unknown. Instead,
this just gives a common syntax for receiving this pair with some default
values.
One can then construct a `Pin` using the `@named` helper macro:

```@example acausal
@named mypin1 = Pin()
```

Next, we build our ground node. A ground node is just a pin that is connected
to a constant voltage reservoir, typically taken to be `V = 0`. Thus to define
this component, we generate a `System` with a `Pin` subcomponent and specify
that the voltage in such a `Pin` is equal to zero. This gives:

```@example acausal
function Ground(; name)
    @named g = Pin()
    eqs = [g.v ~ 0]
    System(eqs, t, [], []; systems = [g], name)
end
```

Next we build a `OnePort`: an abstraction for all simple electric component
with two pins. The voltage difference between the positive pin and the negative
pin is the voltage of the component, the current between two pins must sum to
zero, and the current of the component equals to the current of the positive
pin.

```@example acausal
function OnePort(; name)
    @named p = Pin()
    @named n = Pin()
    @variables v(t) i(t)
    eqs = [
        v ~ p.v - n.v
        0 ~ p.i + n.i
        i ~ p.i
    ]
    System(eqs, t, [v, i], []; systems = [p, n], name)
end
```

Next we build a resistor. A resistor is an object that has two `Pin`s, the positive
and the negative pins, and follows Ohm's law: `v = i*r`. The voltage of the
resistor is given as the voltage difference across the two pins, while by conservation
of charge we know that the current in must equal the current out, which means
(no matter the direction of the current flow) the sum of the currents must be
zero. This leads to our resistor equations:

```@example acausal
function Resistor(; name, R = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    @parameters R = R # Sets the default resistance
    eqs = [v ~ i * R]
    extend(System(eqs, t, [], [R]; name), oneport)
end
```

Notice that we have created this system with a default parameter `R` for the
resistor's resistance. By doing so, if the resistance of this resistor is not
overridden by a higher level default or overridden at `ODEProblem` construction
time, this will be the value of the resistance. Also, note the use of `extend`.
For the `Resistor`, we want to simply inherit `OnePort`'s
equations and unknowns and extend them with a new equation. Note that `v`, `i` are unpacked
from `oneport` using `@unpack`.

Using our knowledge of circuits, we similarly construct the `Capacitor`:

```@example acausal
function Capacitor(; name, C = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    @parameters C = C
    eqs = [D(v) ~ i / C]
    extend(System(eqs, t, [], [C]; name), oneport)
end
```

Now we want to build a constant voltage electric source term. We can think of
this as similarly being a two pin object, where the object itself is kept at a
constant voltage, essentially generating the electric current. We would then
model this as:

```@example acausal
function ConstantVoltage(; name, V = 1.0)
    @named oneport = OnePort()
    @unpack v = oneport
    @parameters V = V
    eqs = [V ~ v]
    extend(System(eqs, t, [], [V]; name), oneport)
end
```

### Connecting and Simulating Our Electric Circuit

Now we are ready to simulate our circuit. Let's build our four components:
a `resistor`, `capacitor`, `source`, and `ground` term. For simplicity, we will
make all of our parameter values 1.0, except the resistor which we set to 2.0.

```@example acausal
@named resistor = Resistor(R = 2.0)
@named capacitor = Capacitor(C = 1.0)
@named source = ConstantVoltage(V = 1.0)
@named ground = Ground()

rc_eqs = [
    connect(source.p, resistor.p)
    connect(resistor.n, capacitor.p)
    connect(capacitor.n, source.n)
    connect(capacitor.n, ground.g)
]

@named rc_model = System(rc_eqs, t; systems = [resistor, capacitor, source, ground])
rc_model = mtkcompile(rc_model)
```

This model is acausal because we have not specified anything about the causality of the model. We have
simply specified what is true about each of the variables. This forms a system
of differential-algebraic equations (DAEs) which define the evolution of each
unknown of the system. The equations are:

```@example acausal
equations(expand_connections(rc_model))
```

the unknowns are:

```@example acausal
unknowns(rc_model)
```

and the parameters are:

```@example acausal
parameters(rc_model)
```

The observed equations are:

```@example acausal
observed(rc_model)
```

## Solving this System

We are left with a system of only two equations
with two unknown variables. One of the equations is a differential equation,
while the other is an algebraic equation. We can then give the values for the
initial conditions of our unknowns, and solve the system by converting it to
an ODEProblem in mass matrix form and solving it with an [ODEProblem mass matrix
DAE solver](https://docs.sciml.ai/DiffEqDocs/stable/solvers/dae_solve/#OrdinaryDiffEq.jl-(Mass-Matrix)).
This is done as follows:

```@example acausal
u0 = [rc_model.capacitor.v => 0.0]

prob = ODEProblem(rc_model, u0, (0, 10.0))
sol = solve(prob)
plot(sol)
```

By default, this plots only the unknown variables that had to be solved for.
However, what if we wanted to plot the timeseries of a different variable? Do
not worry, that information was not thrown away! Instead, transformations
like `mtkcompile` simply change unknown variables into observables which are
defined by `observed` equations.

```@example acausal
observed(rc_model)
```

These are explicit algebraic equations which can then be used to reconstruct
the required variables on the fly. This leads to dramatic computational savings
because implicitly solving an ODE scales like O(n^3), so making there be as
few unknowns as possible is good!

The solution object can be accessed via its symbols. For example, let's retrieve
the voltage of the resistor over time:

```@example acausal
sol[rc_model.resistor.v]
```

or we can plot the timeseries of the resistor's voltage:

```@example acausal
plot(sol, idxs = [rc_model.resistor.v])
```

Although it may be more confusing than helpful here, we can of course also plot all unknown and observed variables together:

```@example acausal
plot(sol, idxs = [unknowns(rc_model); observables(rc_model)])
```
