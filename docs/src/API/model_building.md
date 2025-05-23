# Model building reference

This page lists functionality and utilities related to building hierarchical models. It is
recommended to read the page on the [`System`](@ref System_type) before this.

## Hierarchical model composition

The `System` data structure can represent a tree-like hierarchy of systems for building models
from composable blocks. The [`ModelingToolkit.get_systems`](@ref) function can be used for
querying the subsystems of a system. The `@component` macro should be used when writing
building blocks for model composition.

```@docs
@component
```

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

### Connection semantics

ModelingToolkit implements connection semantics similar to those in the [Modelica specification](https://specification.modelica.org/maint/3.6/connectors-and-connections.html).
We do not support the concept of `inner` and `outer` elements or `expandable` connectors.
Connectors in ModelingToolkit are systems with the appropriate metadata added via the `@connector`
macro.

```@docs
Symbolics.connect
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
Equality
Stream
Flow
```

These are specified using the `connect` metadata. ModelingToolkit also supports `instream`

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
