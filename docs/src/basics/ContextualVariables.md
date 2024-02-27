# Contextual Variable Types

ModelingToolkit.jl has a system of contextual variable types which allows for
helping the system transformation machinery do complex manipulations and
automatic detection. The standard variable definition in ModelingToolkit.jl is
the `@variable` which is defined by
[Symbolics.jl](https://docs.sciml.ai/Symbolics/stable/). For example:

```julia
@variables x y(x)
```

This is used for the “normal” variable of a given system, like the unknowns of a
differential equation or objective function. All the macros below support
the same syntax as `@variables`.

## Parameters

All modeling projects have some form of parameters. `@parameters` marks a variable
as being the parameter of some system, which allows automatic detection algorithms
to ignore such variables when attempting to find the unknowns of a system.

## [Constants](@id constants)

Constants, defined by e.g. `@constants myconst1` are like parameters that:

  - always have a default value, which must be assigned when the constants are
    declared
  - do not show up in the list of parameters of a system.

The intended use-cases for constants are:

  - representing literals (e.g., π) symbolically, which results in cleaner
    Latexification of equations (avoids turning `d ~ 2π*r` into `d = 6.283185307179586 r`)
  - allowing auto-generated unit conversion factors to live outside the list of
    parameters
  - representing fundamental constants (e.g., speed of light `c`) that should never
    be adjusted inadvertently.

## Wildcard Variable Arguments

```julia
@variables u(..)
```

It is possible to define a dependent variable which is an open function as above,
for which its arguments must be specified each time it is used. This is useful with
PDEs for example, where one may need to use `u(t, x)` in the equations, but will
need to be able to write `u(t, 0.0)` to define a boundary condition at `x = 0`.

## Variable metadata

In many engineering systems, some variables act like “flows” while others do not.
For example, in circuit models you have current which flows, and the related
voltage which does not. Or in thermal models you have heat flows. In these cases,
the `connect` statement enforces conservation of flow between all of the connected
components.

For example, the following specifies that `x` is a 2x2 matrix of flow variables
with the unit m^3/s:

```julia
@variables x[1:2, 1:2] [connect = Flow; unit = u"m^3/s"]
```

ModelingToolkit defines `connect`, `unit`, `noise`, and `description` keys for
the metadata. One can get and set metadata by

```julia
julia> @variables x [unit = u"m^3/s"];

julia> hasmetadata(x, VariableUnit)
true

julia> ModelingToolkit.get_unit(x)
m³ s⁻¹

julia> x = setmetadata(x, VariableUnit, u"m/s")
x

julia> ModelingToolkit.get_unit(x)
m s⁻¹
```

See [Symbolic Metadata](@ref symbolic_metadata) for more details on variable metadata.
