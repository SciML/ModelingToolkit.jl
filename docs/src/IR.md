# ModelingToolkit IR

ModelingToolkit IR, which falls under the `Expression` abstract type, mirrors
the Julia AST but allows for easy mathematical manipulation by itself following
mathematical semantics. The base of the IR is the `Variable` type, which defines
a symbolic variable. These variables are combined using `Operation`s, which are
registered functions applied to the various variables. These `Operation`s then
perform automatic tracing, so normal mathematical functions applied to an `Operation`
generate a new `Operation`. For example, `op1 = x+y` is one `Operation` and
`op2 = 2z` is another, and so `op1*op2` is another `Operation`. Then, at the top,
an `Equation`, normally written as `op1 ~ op2`, defines the symbolic equality
between two operations.

### Types

```@docs
Expression
Variable
ModelingToolkit.Constant
Operation
Equation
```

### Function Registration

The ModelingToolkit graph only allowed for registered Julia functions for the
operations. All other functions are automatically traced down to registered
functions. By default, ModelingToolkit.jl pre-registers the common functions
utilized in the AD package ruleset [DiffRules.jl](https://github.com/JuliaDiff/DiffRules.jl)
and pre-defines their derivatives. However, the user can utilize the `@register`
macro to add their function to allowed functions of the computation graph.

```@docs
@register
```

### Derivatives and Differentials

A `Differential(op)` is a partial derivative with respect to the operation `op`,
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
```

Note that the generation of sparse matrices simply follows from the Julia semantics
imbued on the IR, so `sparse(jac)` changes a dense Jacobian to a sparse Jacobian
matrix.

### Adding Derivatives

There is a large amount of derivatives pre-defined by
[DiffRules.jl](https://github.com/JuliaDiff/DiffRules.jl). Note that `Expression`
types are defined as `<:Real`, and thus any functions which allow the use of real
numbers can automatically be traced by the derivative mechanism. Thus, for example:

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
You should return an `Operation` for the derivative of your function.

For example, `sin(t)`'s derivative (by `t`) is given by the following:

```julia
ModelingToolkit.derivative(::typeof(sin), args::NTuple{1,Any}, ::Val{1}) = cos(args[1])
```

### IR Manipulation

ModelingToolkit.jl provides functionality for easily manipulating `Expression`
types. Most of the functionality comes by the `Expression` type obeying the
standard mathematical semantics. For example, if one has `A` a matrix of
`Expression`, then `A^2` calculates the `Expression`s for the squared matrix.
In that sense, it is encouraged that one uses standard Julia for performing a
lot of the manipulation on the IR, as, for example, calculating the sparse form
of the matrix via `sparse(A)` is valid, legible, and easily understandable
to all Julia programmers.

Other additional manipulation functions are given below.

```@docs
simplify_constants
rename
get_variables
substitute_expr!
```
