# Contextual Variable Types

ModelingToolkit.jl has a system of contextual variable types which allows for
helping the system transformation machinery do complex manipulations and
automatic detection. The standard variable definition in ModelingToolkit.jl is
the `@variable` which is defined by
[Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl). For example:

```julia
@variables x y(x)
```

This is used for the "normal" variable of a given system, like the states of a
differential equation or objective function. All of the macros below support
the same syntax as `@variables`.

## Parameters

All modeling projects have some form of parameters. `@parameters` marks a variable
as being the parameter of some system, which allows automatic detection algorithms
to ignore such variables when attempting to find the states of a system.

## Flow Variables (TODO)

In many engineering systems some variables act like "flows" while others do not.
For example, in circuit models you have current which flows, and the related
voltage which does not. Or in thermal models you have heat flows. In these cases,
the `connect` statement enforces conservation of flow between all of the connected
components.

For example, the following specifies that `x` is a 2x2 matrix of flow variables:

```julia
@variables x[1:2,1:2] type=flow
```

## Stream Variables

TODO

## Brownian Variables

TODO
