# Callable parameters and interpolating data

ModelingToolkit.jl allows creating parameters that represent functions to be called. This
is especially useful for including interpolants and/or lookup tables inside ODEs. In this
tutorial we will create an `ODESystem` which employs callable parameters to interpolate data
inside an ODE and go over the various syntax options and their implications.

## Callable parameter syntax

The syntax for callable parameters declared via `@parameters` must be one of the following

 1. `(fn::FType)(..)`
 2. `fn(::argType1, ::argType2, ...)`

In the first case, the parameter is callable with any number/combination of arguments, and
has a type of `FType` (the callable must be a subtype of `FType`). In the second case,
the parameter is callable with as many arguments as declared, and all must match the
declared types.

By default, the return type of the callable symbolic is inferred to be `Real`. To change
this, a `::retType` annotation can be added at the end.

To declare a function that returns an array of values, the same array syntax can be used
as for normal variables:

```julia
@parameters (foo::FType)(..)[1:3]::retType
@parameters foo(::argType1, ::argType2)[1:3]::retType
```

`retType` here is the `eltype` of the returned array.

## Storage of callable parameters

Callable parameters declared with the `::FType` syntax will be stored in a `Vector{FType}`.
Thus, if `FType` is non-concrete, the buffer will also be non-concrete. This is sometimes
necessary to allow the value of the callable to be switched out for a different type without
rebuilding the model. Typically this syntax is preferable when `FType` is concrete or
a small union.

Callable parameters declared with the `::argType1, ...` syntax will be stored in a
`Vector{FunctionWrappers.FunctionWrapper{retType, Tuple{argType1, ...}}}`. This suffers
the small overhead of a `FunctionWrapper` and restricts the signature of the callable,
symbolic, but allows storing the parameter in a type-stable manner and swapping it out.
This is preferable when the values that the callable can take do not share a common
subtype. For example, when a callable can represent the activation of a neural network
and can be `tanh`, `sigmoid`, etc. which have a common ancestor of `Function`.

If both `::FType` and `::argType`s are specified, `::FType` takes priority. For example,

```julia
@parameters (p::LinearInterpolation)(::Real)
```

`p` will be stored in a `Vector{LinearInterpolation}`. If `::LinearInterpolation` was not
specified, it would be stored in a `Vector{FunctionWrapper{Real, Tuple{Real}}}`.

## Example using interpolations

```@example callable
using ModelingToolkit
using OrdinaryDiffEq
using DataInterpolations
using ModelingToolkit: t_nounits as t, D_nounits as D

ts = collect(0.0:0.1:10.0)
spline = LinearInterpolation(ts .^ 2, ts)
Tspline = typeof(spline)
@variables x(t)
@parameters (interp::Tspline)(..)

@mtkbuild sys = ODESystem(D(x) ~ interp(t), t)
```

The derivative of `x` is obtained via an interpolation from DataInterpolations.jl. Note
the parameter syntax. The `(..)` marks the parameter as callable. `(interp::Tspline)`
indicates that the parameter is of type `Tspline`.

```@example callable
prob = ODEProblem(sys, [x => 0.0], (0.0, 1.0), [interp => spline])
solve(prob)
```

Note that the the following will not work:

```julia
ODEProblem(
    sys; [x => 0.0], (0.0, 1.0), [interp => LinearInterpolation(0.0:0.1:1.0, 0.0:0.1:1.0)])
```

Since the type of the spline doesn't match.
