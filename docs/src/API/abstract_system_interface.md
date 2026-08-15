# [The `AbstractSystem` Interface](@id abstract_system_interface)

[`ModelingToolkit.AbstractSystem`](@ref) is the supertype of every ModelingToolkit system.
This page states what a value is guaranteed to support purely by virtue of being an
`AbstractSystem`, so that code which inspects, traverses, or reports on systems can be
written against the abstraction instead of against the fields of one concrete type.

The contracts below are written as dispatch signatures with the semantics that must hold.
Where a function is only available when the system stores a particular piece of data, that
is stated explicitly along with the query used to detect availability; everything else
holds for every subtype.

!!! warning "Interface scope"

    This interface covers structural inspection of a system: its name, its subsystem
    hierarchy, its independent variables, and the symbolic collections it stores. The
    numerical pipeline — `mtkcompile`, problem and function construction, code generation,
    initialization, and the indexing surface those produce — is defined on the concrete
    [`System`](@ref System_type) type and is not part of the `AbstractSystem` contract. See
    [Outside the interface](@ref abstract_system_interface_boundary).

## Required Implementation

A subtype must supply two pieces of information. Both have a default implementation that
reads a field of the same name, so a subtype which declares those fields implements the
requirement without writing any methods.

| Signature                                                   | Default                | Meaning                            |
|:----------------------------------------------------------- |:---------------------- |:---------------------------------- |
| `Base.nameof(sys::AbstractSystem)::Symbol`                   | reads the `name` field | The name of this system.           |
| `ModelingToolkit.get_systems(sys::AbstractSystem)::Vector`   | reads the `systems` field | The direct subsystems of this system. |

The following rules apply to every subtype:

  - [`nameof`](@ref) returns a `Symbol` and must be unique among the systems returned by
    [`ModelingToolkit.get_systems`](@ref) of a common parent. Namespacing uses this name as
    the prefix, so duplicate sibling names produce colliding variable names rather than an
    error.
  - [`ModelingToolkit.get_systems`](@ref) returns only the *direct* children, each of which
    is itself an `AbstractSystem`. It never returns transitive descendants.
  - The hierarchy reachable through [`ModelingToolkit.get_systems`](@ref) must be finite and
    acyclic. Every generic accessor below recurses through it without cycle detection.
  - Systems are treated as immutable values. A generic function may retain, compare, or hash
    a system, and may assume that a system it received earlier still describes the same
    model.

## Optional Storage and the `has_x`/`get_x` Protocol

Everything beyond the name and the subsystem list is optional storage. ModelingToolkit
recognizes a fixed set of field names — those listed among
[the accessor functions of `System`](@ref System_type) — and defines an accessor pair for
each:

```julia
ModelingToolkit.has_x(sys::AbstractSystem)::Bool
ModelingToolkit.get_x(sys::AbstractSystem)
```

  - `has_x` is defined for every `AbstractSystem` and never throws. It reports whether this
    system stores `x` at all; it says nothing about whether the stored value is empty.
  - `get_x` returns the value stored *locally* by `sys`, without namespacing and without
    consulting subsystems. It throws if the system does not store `x`, so generic code must
    guard it with `has_x`. The event accessors
    [`ModelingToolkit.get_continuous_events`](@ref) and
    [`ModelingToolkit.get_discrete_events`](@ref) are the documented exceptions: they return
    an empty vector instead of throwing.
  - The pair exists only for recognized field names. A subtype is free to hold additional
    data under other field names, but ModelingToolkit will not expose it, and those fields
    take no part in this interface.

Because optional storage is genuinely optional, a generic function that needs a piece of
data must decide what to do when it is absent, rather than assuming the layout of
[`System`](@ref System_type).

## Accessors Available for Every Subtype

These functions are total on `AbstractSystem`: they are defined for every subtype and
substitute a documented value when the storage they would read is absent.

| Signature                                                                       | Value when the backing storage is absent |
|:------------------------------------------------------------------------------- |:---------------------------------------- |
| `ModelingToolkit.independent_variable(sys::AbstractSystem)`                       | `nothing`                                |
| `ModelingToolkit.independent_variables(sys::AbstractSystem)::Vector`              | empty vector                             |
| `ModelingToolkit.description(sys::AbstractSystem)::String`                        | `""`                                     |
| `ModelingToolkit.iscomplete(sys::AbstractSystem)::Bool`                           | `false`                                  |
| `ModelingToolkit.does_namespacing(sys::AbstractSystem)::Bool`                     | `!iscomplete(sys)`                       |
| `ModelingToolkit.assertions(sys::AbstractSystem)::Dict`                           | empty dictionary                         |
| `continuous_events(sys::AbstractSystem)::Vector`                                  | empty vector                             |
| `discrete_events(sys::AbstractSystem)::Vector`                                    | empty vector                             |

[`ModelingToolkit.independent_variable`](@ref ModelingToolkitBase.independent_variable) is
the scalar extension point: a subtype that stores its independent variable under a different
field name, or computes it, adds a method to that function. The vector-valued
[`ModelingToolkit.independent_variables`](@ref ModelingToolkitBase.independent_variables) is
derived from it and must not be extended.

ModelingToolkit also implements the parts of the
[SymbolicIndexingInterface](https://docs.sciml.ai/SymbolicIndexingInterface/stable/) that
follow from the above for every subtype:

```julia
SymbolicIndexingInterface.independent_variable_symbols(sys::AbstractSystem)::Vector
SymbolicIndexingInterface.constant_structure(sys::AbstractSystem)::Bool
```

## Accessors Backed by Optional Storage

Each function below reads one optional field and therefore requires that the system stores
it. Guard the call with the corresponding `has_x` query.

| Accessor                                                                                                                 | Required storage       |
|:------------------------------------------------------------------------------------------------------------------------ |:---------------------- |
| [`equations`](@ref), [`ModelingToolkit.equations_toplevel`](@ref), [`alg_equations`](@ref), [`diff_equations`](@ref), [`has_alg_equations`](@ref), [`has_diff_equations`](@ref), [`ModelingToolkit.namespace_equations`](@ref) | `eqs`                  |
| [`unknowns`](@ref), [`ModelingToolkit.unknowns_toplevel`](@ref)                                                            | `unknowns`             |
| [`parameters`](@ref), [`ModelingToolkit.parameters_toplevel`](@ref)                                                        | `ps`                   |
| [`observed`](@ref), [`observables`](@ref)                                                                                  | `observed`             |
| [`initial_conditions`](@ref)                                                                                               | `initial_conditions`   |
| [`guesses`](@ref)                                                                                                          | `guesses`              |
| [`initialization_equations`](@ref)                                                                                         | `initialization_eqs`   |
| [`constraints`](@ref)                                                                                                      | `constraints`          |
| [`cost`](@ref)                                                                                                             | `costs`, `consolidate` |
| [`brownians`](@ref)                                                                                                        | `brownians`            |
| [`jumps`](@ref)                                                                                                            | `jumps`                |

The variable and parameter queries of the SymbolicIndexingInterface are answered from
[`unknowns`](@ref) and [`parameters`](@ref), so they carry the same storage requirement:

```julia
SymbolicIndexingInterface.variable_symbols(sys::AbstractSystem)
SymbolicIndexingInterface.is_variable(sys::AbstractSystem, sym)
SymbolicIndexingInterface.variable_index(sys::AbstractSystem, sym)
SymbolicIndexingInterface.parameter_symbols(sys::AbstractSystem)
SymbolicIndexingInterface.is_parameter(sys::AbstractSystem, sym)
SymbolicIndexingInterface.parameter_index(sys::AbstractSystem, sym)
```

## Hierarchy and Namespacing Semantics

The accessors in the previous section come in two families, and the difference between them
is the whole of the hierarchy contract.

A **hierarchical** accessor returns the locally stored value followed by the result of the
same accessor on each system in [`ModelingToolkit.get_systems`](@ref), applied recursively
and namespaced at each level. [`equations`](@ref), [`unknowns`](@ref),
[`parameters`](@ref), [`observed`](@ref), [`observables`](@ref),
[`initial_conditions`](@ref), [`guesses`](@ref), [`initialization_equations`](@ref),
[`constraints`](@ref), [`cost`](@ref), [`brownians`](@ref), [`jumps`](@ref),
[`continuous_events`](@ref), [`discrete_events`](@ref), and
[`ModelingToolkit.assertions`](@ref) all behave this way, as do
[`alg_equations`](@ref), [`diff_equations`](@ref), [`has_alg_equations`](@ref) and
[`has_diff_equations`](@ref), which classify the result of [`equations`](@ref):

  - Local entries come first, then subsystems in the order
    [`ModelingToolkit.get_systems`](@ref) returns them. [`parameters`](@ref) additionally
    removes duplicates, since the same parameter may be reachable through more than one
    subsystem.
  - A symbol `x` owned by subsystem `sub` appears as
    [`ModelingToolkit.renamespace`](@ref)`(sub, x)`, which prefixes `nameof(sub)`. Nesting
    composes, so a symbol two levels down carries both prefixes.
  - The result is a fresh collection whenever the system has subsystems. When it has none,
    the collection may alias the system's own storage. Treat every returned collection as
    read-only.

A **toplevel** accessor — [`ModelingToolkit.equations_toplevel`](@ref),
[`ModelingToolkit.unknowns_toplevel`](@ref),
[`ModelingToolkit.parameters_toplevel`](@ref),
[`ModelingToolkit.continuous_events_toplevel`](@ref),
[`ModelingToolkit.discrete_events_toplevel`](@ref) — returns only what the system stores
itself, with no namespacing and no recursion. `ModelingToolkit.get_x` is the same view with
no interpretation applied at all.

## [Outside the interface](@id abstract_system_interface_boundary)

The following are properties of the concrete [`System`](@ref System_type) type rather than of
`AbstractSystem`, and a generic function must not assume them:

  - **Compilation and numerics.** `complete`, `mtkcompile`, every problem and function
    constructor, code generation, and initialization require a `System`.
  - **Display, reconstruction, and metadata.** Printing a system, `Symbolics.rename`,
    `ConstructionBase.setproperties`, `ModelingToolkit.toggle_namespacing`, and the
    `SymbolicUtils` metadata accessors reconstruct the system from its fields and depend on
    the internal representation used by `System`.
  - **`getproperty` for subsystems.** `sys.subsystem` re-namespaces the retrieved subsystem
    and therefore reconstructs it. Reach subsystems generically through
    [`ModelingToolkit.get_systems`](@ref) instead.
  - **The remaining SymbolicIndexingInterface surface.** `is_time_dependent`,
    `all_symbols`, `all_variable_symbols`, `default_values`, and the observed-function
    machinery are implemented for `System`.

A package that needs one of these should build a [`System`](@ref System_type) rather than
introduce a subtype, and should route data through the documented `System` constructors.

## Example

The system below implements the required part of the interface by declaring the `name` and
`systems` fields, and opts into the equation, unknown, and parameter accessors by declaring
those fields too.

```@example abstract_system_interface
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

struct ComponentSystem <: ModelingToolkit.AbstractSystem
    eqs::Vector{Equation}
    unknowns::Vector
    ps::Vector
    iv::Any
    name::Symbol
    systems::Vector
end

@variables x(t) y(t)
@parameters p q

leaf = ComponentSystem([D(x) ~ p * x], [x], [p], t, :leaf, [])
root = ComponentSystem([D(y) ~ q * y], [y], [q], t, :root, [leaf])

equations(root)
```

The hierarchical accessors namespace the contents of `leaf` under its name, while the
toplevel accessors report only what `root` stores:

```@example abstract_system_interface
unknowns(root), ModelingToolkit.unknowns_toplevel(root)
```

Storage the subtype does not declare is reported as absent rather than as empty, and the
accessor it backs must not be called. Here `iv` is declared and `observed` is not:

```@example abstract_system_interface
ModelingToolkit.has_iv(root), ModelingToolkit.has_observed(root)
```

## Validating a Subtype

Exercise the subtype through the functions on this page only. A test that reads a field
directly, or that calls an entry point listed under
[Outside the interface](@ref abstract_system_interface_boundary), is testing the `System`
implementation rather than the interface.

`lib/ModelingToolkitBase/test/abstractsystem_interface.jl` is the executable form of this
page: it defines minimal subtypes and asserts each rule stated here.
