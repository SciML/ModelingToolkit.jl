# ModelingToolkit IR

ModelingToolkit IR mirrors the Julia AST but allows for easy mathematical
manipulation by itself following mathematical semantics. The base of the IR is
the `Sym` type, which defines a symbolic variable. Registered (mathematical)
functions on `Sym`s (or `Term`s) return `Term`s.  For example, `op1 = x+y` is
one `Term` and `op2 = 2z` is another, and so `op1*op2` is another `Term`. Then,
at the top, an `Equation`, normally written as `op1 ~ op2`, defines the
symbolic equality between two operations.

### Types

```@docs
Sym
Term
Equation
```

### A note about functions restricted to `Number`s

`Sym` and `Term` objects are NOT subtypes of `Number`. ModelingToolkit provides
a simple wrapper type called `Num` which is a subtype of `Real`. `Num` wraps
either a Sym or a Term or any other object, defines the same set of operations
as symbolic expressions and forwards those to the values it wraps. You can use
`ModelingToolkit.value` function to unwrap a `Num`.

By default, the `@variables` and `@parameters` functions return Num-wrapped
objects so as to allow calling functions which are restricted to `Number` or
`Real`.

### Function Registration

The ModelingToolkit graph only allowed for registered Julia functions for the
operations. All other functions are automatically traced down to registered
functions. By default, ModelingToolkit.jl pre-registers the common functions
utilized in [SymbolicUtils.jl](https://github.com/JuliaSymbolics/SymbolicUtils.jl)
and pre-defines their derivatives. However, the user can utilize the `@register`
macro to add their function to allowed functions of the computation graph.

```@docs
@register
```

### Derivatives and Differentials

A `Differential(op)` is a partial derivative with respect to `op`,
which can then be applied to some other operations. For example, `D=Differential(t)`
is what would commonly be referred to as `d/dt`, which can then be applied to
other operations using its function call, so `D(x+y)` is `d(x+y)/dt`.

By default, the derivatives are left unexpanded to capture the symbolic
representation of the differential equation. If the user would like to expand
out all of the differentials, the `expand_derivatives` function eliminates all
of the differentials down to basic one-variable expressions.

```@docs
ModelingToolkit.derivative
Differential
expand_derivatives
ModelingToolkit.jacobian
ModelingToolkit.gradient
ModelingToolkit.hessian
```

For jacobians which are sparse, use the `sparsejacobian` function.
For hessians which are sparse, use the `sparsehessian` function.

### Adding Derivatives

There is a large amount of derivatives pre-defined by
[DiffRules.jl](https://github.com/JuliaDiff/DiffRules.jl).

```julia
f(x,y,z) = x^2 + sin(x+y) - z
```

automatically has the derivatives defined via the tracing mechanism. It will do
this by directly building the operation the internals of your function and
differentiating that.

However, in many cases you may want to define your own derivatives so that way
automatic Jacobian etc. calculations can utilize this information. This can
allow for more succinct versions of the derivatives to be calculated in order
to better scale to larger systems. You can define derivatives for your own
function via the dispatch:

```julia
# `N` arguments are accepted by the relevant method of `my_function`
ModelingToolkit.derivative(::typeof(my_function), args::NTuple{N,Any}, ::Val{i})
```

where `i` means that it's the derivative with respect to the `i`th argument. `args` is the
array of arguments, so, for example, if your function is `f(x,t)`, then `args = [x,t]`.
You should return an `Term` for the derivative of your function.

For example, `sin(t)`'s derivative (by `t`) is given by the following:

```julia
ModelingToolkit.derivative(::typeof(sin), args::NTuple{1,Any}, ::Val{1}) = cos(args[1])
```

### IR Manipulation

ModelingToolkit.jl provides functionality for easily manipulating expressions.
Most of the functionality comes by the expression objects obeying the standard
mathematical semantics. For example, if one has `A` a matrix of symbolic
expressions wrapped in `Num`, then `A^2` calculates the expressions for the
squared matrix.  In that sense, it is encouraged that one uses standard Julia
for performing a lot of the manipulation on the IR, as, for example,
calculating the sparse form of the matrix via `sparse(A)` is valid, legible,
and easily understandable to all Julia programmers.

Other additional manipulation functions are given below.

```@docs
simplify_constants
rename
get_variables
substitute_expr!
```
