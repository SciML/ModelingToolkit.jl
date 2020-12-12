# ModelingToolkit IR

ModelingToolkit IR mirrors the Julia AST but allows for easy mathematical
manipulation by itself following mathematical semantics.

The IR is made of the following types:

1. `Sym`:
  a. `Sym{T}(:x)` (created using `@variables x::T`) represents a variable of type `T`. If `::T` is omitted, defaults to `Real`
  b. `Sym{Parameter{T}}(:ρ)` (created using `@variables ρ::T`) represents a _parameter_ of type `T`.
  c. `Sym{FnType{Tuple{X, Y}, Z}}(:f)` (created using `@variables f(::X, ::Y)::Z`) represents a variable which behaves as a function which takes 2 arguments of symbolic type X and Y respectively and returns an object of symbolic type `Z`. 
  Supports: [`nameof`](@ref), [`symtype`](@ref)
2. `Term`:
  a. when a mathematical operation is called on a `Sym` of the first two kind above (variable and parameter),it results in a `Term`. `Term`s are also closed under the same mathematical operations.
  b. when a symbolic function (an object of the 3rd kind of Sym decsribed above) is called with arguments of the appropriate type, it causes a `Term` to be created with the `Sym` as its `operation`.
  Supports: [`operation`](@ref), [`arguments`](@ref)
3. `Symbolic`: the super type of `Sym` and `Term`, used for convenience.
4. `Num`: wraps either a `Symbolic` or a `Real` and is itself a subtype of `Real`. When mathematical operations are called with one or more `Num`s the results are computed with the unwrapped values and then wrapped with `Num`. The purpose of `Num` is to allow ModelingToolkit expressions to propagate through code that is restricted to `Real`s and also to make it possible for correct handling of `one` and `zero` required by common operations on Arrays of expressions. (e.g. `sum(Num[]) == 0`; `qr(num_matrix)` etc.).
  Supports: [`value`](@ref)

`Expression` is a type alias for `Union{Term{<:Real}, Sym{<:Real}, Num, Real}`. This refers to any numerically typed symbolic or literal value. This is mainly used for the purpose of documentation.

User facing APIs (e.g. `jacobian`, `@variables`) all take and return `Num`s, while most of ModelingToolkit internals work on the unwrapped expressions, namely `Sym` and `Term`s.

```@docs
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

```julia
julia> @parameters t; @variables x y z(t);

julia> ModelingToolkit.operation(ModelingToolkit.value(x + y))
+ (generic function with 377 methods)

julia> ModelingToolkit.operation(ModelingToolkit.value(z))
z(::Any)::Real

julia> ModelingToolkit.arguments(ModelingToolkit.value(x + y))
2-element Vector{Sym{Real}}:
 x
 y
```

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
get_variables
substitute
toparam
tosymbol
makesym
diff2term
```
