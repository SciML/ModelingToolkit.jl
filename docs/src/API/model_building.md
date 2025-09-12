# [Model building reference](@id model_building_api)

This page lists functionality and utilities related to building hierarchical models. It is
recommended to read the page on the [`System`](@ref System_type) before this.

## Common definitions of `t` and `D`

ModelingToolkit provides common definitions for the independent variable `t` (time) and the
derivative with respect to it `D`.

```@docs
ModelingToolkit.t_nounits
ModelingToolkit.D_nounits
ModelingToolkit.t
ModelingToolkit.D
ModelingToolkit.t_unitful
ModelingToolkit.D_unitful
```

Users are recommended to use the appropriate common definition in their models. The required
definitions can be imported with convenient aliased names. For example:

```julia
using ModelingToolkit: t_nounits as t, D_nounits as D
```

Allows using `t` and `D` to refer to `t_nounits` and `D_nounits` respectively.

## Hierarchical model composition

The `System` data structure can represent a tree-like hierarchy of systems for building models
from composable blocks. The [`ModelingToolkit.get_systems`](@ref) function can be used for
querying the subsystems of a system. The `@component` macro should be used when writing
building blocks for model composition.

```@docs
@component
```

Every constructor function should build either a component or a connector. Components define
the dynamics of the system. Connectors are used to connect components together and propagate
information between them. See also [`@connector`](@ref).

### Scoping of variables

When building hierarchical systems, is is often necessary to pass variables from a parent system
to the subsystems. If done naively, this will result in the child system assuming it "owns" the
variables passed to it and any occurrences of those variables in the child system will be
namespaced. To prevent this, ModelingToolkit has the concept of variable scope. The scope allows
specifying which system a variable belongs to relative to the system in which it is used.

```@docs
LocalScope
ParentScope
GlobalScope
```

Note that the scopes must be applied to _individual variables_ and not expressions. For example,
`ParentScope(x + y)` is incorrect. Instead, `ParentScope(x) + ParentScope(y)` is the correct usage.
Applying the same scope (more generally, the same function) to all variables in an expression is a
common task, and ModelingToolkit exposes a utility for the same:

```@docs
ModelingToolkit.apply_to_variables
```

It is still tedious to manually use `apply_to_variables` on any symbolic expression passed to a
subsystem. The `@named` macro automatically wraps all symbolic arguments in `ParentScope` and
uses the identifier being assigned as the name of the system.

```@docs
@named
```

### Exploring the tree structure

The `System` type implements the `AbstractTrees` interface. This can be used to explore the
hierarchical structure.

```@docs
hierarchy
```

### [Connection semantics](@id connect_semantics)

ModelingToolkit implements connection semantics similar to those in the [Modelica specification](https://specification.modelica.org/maint/3.6/connectors-and-connections.html).
We do not support the concept of `inner` and `outer` elements or `expandable` connectors.
Connectors in ModelingToolkit are systems with the appropriate metadata added via the `@connector`
macro.

```@docs
connect
domain_connect
@connector
```

Connections can be expanded using `expand_connections`.

```@docs
expand_connections
```

Similar to the `stream` and `flow` keyword arguments in the specification, ModelingToolkit
allows specifying how variables in a connector behave in a connection.

```@docs
ModelingToolkit.Equality
Flow
Stream
```

These are specified using the `connect` metadata. ModelingToolkit also supports `instream`.
Refer to the Modelica specification on [Stream connectors](https://specification.modelica.org/maint/3.6/stream-connectors.html)
for more information.

```@docs
instream
```

### System composition utilities

```@docs
extend
compose
substitute_component
```

### Flattening systems

The hierarchical structure can be flattened. This operation is performed during simplification.

```@docs
flatten
```

## System simplification

`System`s can be simplified to reformulate them in a way that enables it to be solved numerically,
and also perform other optimizations. This is done via the `mtkcompile` function. Connection expansion
and flattening are preprocessing steps of simplification.

```@docs
mtkcompile
@mtkcompile
```

It is also possible (though not always advisable) to build numerical problems from systems without
passing them through `mtkcompile`. To do this, the system must first be marked as "complete" via
the `complete` function. This process is used to indicate that a system will not be modified
further and allows ModelingToolkit to perform any necessary preprocessing to it. `mtkcompile`
calls `complete` internally.

```@docs
complete
```

### Exploring the results of simplification

Similar to how [`full_equations`](@ref) returns the equations of a system with all variables
eliminated during `mtkcompile` substituted, we can perform this substitution on an arbitrary
expression.

```@docs
ModelingToolkit.substitute_observed
ModelingToolkit.empty_substitutions
ModelingToolkit.get_substitutions
```

### Experimental simplification

ModelingToolkit may have a variety of experimental simplification passes. These are not
enabled by default, but can be used by passing to the `additional_passes` keyword argument
of `mtkcompile`.

```@docs
ModelingToolkit.IfLifting
```

## Event handling

Time-dependent systems may have several events. These are used to trigger discontinuities
in the model. They compile to standard callbacks from `DiffEqCallbacks.jl`.

```@docs
ModelingToolkit.SymbolicContinuousCallback
ModelingToolkit.SymbolicDiscreteCallback
```

The affect functions for the above callbacks can be symbolic or user-defined functions.
Symbolic affects are handled using equations as described in the [Events](@ref events)
section of the documentation. User-defined functions can be used via `ImperativeAffect`.

```@docs
ModelingToolkit.ImperativeAffect
```

## Modelingtoolkitize

ModelingToolkit can take some numerical problems created non-symbolically and build a
symbolic representation from them.

```@docs
modelingtoolkitize
```

## Using FMUs

ModelingToolkit is capable of importing FMUs as black-box symbolic models. Currently only
a subset of FMU features are supported. This functionality requires importing `FMI.jl`.

```@docs
ModelingToolkit.FMIComponent
```

## Model transformations

ModelingToolkit exposes a variety of transformations that can be applied to models to aid in
symbolic analysis.

```@docs
liouville_transform
fractional_to_ordinary
linear_fractional_to_ordinary
change_of_variables
stochastic_integral_transform
Girsanov_transform
change_independent_variable
add_accumulations
noise_to_brownians
convert_system_indepvar
subset_tunables
respecialize
```

## Hybrid systems

Hybrid systems are dynamical systems involving one or more discrete-time subsystems. These
discrete time systems follow clock semantics - they are synchronous systems and the relevant
variables are only defined at points where the clock ticks.

While ModelingToolkit is unable to simplify, compile and solve such systems on its own, it
has the ability to represent them. Compilation strategies can be implemented independently
on top of [`mtkcompile`](@ref) using the `additional_passes` functionality.

!!! warn
    
    These operators are considered experimental API.

```@docs; canonical = false
Sample
Hold
SampleTime
```

ModelingToolkit uses the clock definition in SciMLBase

```@docs
SciMLBase.TimeDomain
SciMLBase.Clock
SciMLBase.SolverStepClock
SciMLBase.Continuous
```

### State machines

While ModelingToolkit has the capability to represent state machines, it lacks the ability
to compile and simulate them.

!!! warn
    
    This functionality is considered experimental API

```@docs
initial_state
transition
activeState
entry
ticksInState
timeInState
```
