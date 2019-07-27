# ModelingToolkit.jl

[![Build Status](https://travis-ci.org/JuliaDiffEq/ModelingToolkit.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/ModelingToolkit.jl)
[![Coverage Status](https://coveralls.io/repos/JuliaDiffEq/ModelingToolkit.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JuliaDiffEq/ModelingToolkit.jl?branch=master)
[![codecov.io](http://codecov.io/github/JuliaDiffEq/ModelingToolkit.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaDiffEq/ModelingToolkit.jl?branch=master)

ModelingToolkit.jl is an intermediate representation (IR) of computational graphs
for scientific computing problems. Its purpose is to be a common target for
modeling DSLs in order to allow for a common platform for model inspection and
transformation. It uses a tagged variable IR in order to allow specification of
complex models and allow for transformations of models. It has ways to plug into
its function registration and derivative system so that way it can interact
nicely with user-defined routines. Together, this is an abstract form of a
scientific model that is easy for humans to generate but also easy for programs
to manipulate.

## Introduction by Examples

### Example: ODE

Let's build an ODE. First we define some variables. In a differential equation
system, we need to differentiate between our (dependent) variables
and parameters. Therefore we label them as follows:

```julia
using ModelingToolkit

# Define some variables
@parameters t σ ρ β
@variables x(t) y(t) z(t)
@derivatives D'~t
```

Then we build the system:

```julia
eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]
```

Each operation builds an `Operation` type, and thus `eqs` is an array of
`Operation` and `Variable`s. This holds a tree of the full system that can be
analyzed by other programs. We can turn this into a `ODESystem` via:

```julia
de = ODESystem(eqs)
```

where we tell it the variable types and ordering in the first version, or let it
automatically determine the variable types in the second version.
This can then generate the function. For example, we can see the
generated code via:

```julia
generate_function(de, [x,y,z], [σ,ρ,β])

## Which returns:
:((##363, u, p, t)->begin
          let (x, y, z, σ, ρ, β) = (u[1], u[2], u[3], p[1], p[2], p[3])
              ##363[1] = σ * (y - x)
              ##363[2] = x * (ρ - z) - y
              ##363[3] = x * y - β * z
          end
      end)
```

and get the generated function via:

```julia
f = ODEFunction(de, [x,y,z], [σ,ρ,β])
```

### Example: Nonlinear System

We can also build nonlinear systems. Let's say we wanted to solve for the steady
state of the previous ODE. This is the nonlinear system defined by where the
derivatives are zero. We use (unknown) variables for our nonlinear system.

```julia
@variables x y z
@parameters σ ρ β

# Define a nonlinear system
eqs = [0 ~ σ*(y-x),
       0 ~ x*(ρ-z)-y,
       0 ~ x*y - β*z]
ns = NonlinearSystem(eqs, [x,y,z])
nlsys_func = generate_function(ns, [x,y,z], [σ,ρ,β])
```

which generates:

```julia
:((##364, u, p)->begin
          let (x, y, z, σ, ρ, β) = (u[1], u[2], u[3], p[1], p[2], p[3])
              ##364[1] = σ * (y - x)
              ##364[2] = x * (ρ - z) - y
              ##364[3] = x * y - β * z
          end
      end)
```

We can use this to build a nonlinear function for use with NLsolve.jl:

```julia
f = @eval eval(nlsys_func)
# Make a closure over the parameters for for NLsolve.jl
f2 = (du,u) -> f(du,u,(10.0,26.0,2.33))
```

### Example: Arrays of variables

Sometimes it is convenient to define arrays of variables to model things like `x₁,…,x₃`.
The `@variables` and `@parameters` macros support this with the following syntax:

```julia
julia> @variables x[1:3];
julia> x
3-element Array{Operation,1}:
 x₁()
 x₂()
 x₃()

# support for arbitrary ranges and tensors
julia> @variables y[2:3,1:5:6];
julia> y
2×2 Array{Operation,2}:
    y₂̒₁() y₂̒₆()
    y₃̒₁() y₃̒₆()

# also works for dependent variables
julia> @parameters t; @variables z[1:3](t);
julia> z
3-element Array{Operation,1}:
 z₁(t())
 z₂(t())
 z₃(t())
 ```

## Core Principles

The core idea behind ModelingToolkit.jl is that mathematical equations require
context, and thus any symbolic manipulations and full model specifications
requires the ability to handle such context. When writing DSLs, this fact
comes to light very quickly. Every DSL seems to lower to some intermediate
representation from which the final result is computed, but this process means
there's a lot of repeated ideas for every DSL that creates scientific computing
objects like differential equations and nonlinear systems. By having a single
common contexualized IR, this gives DSLs a target to write to so that way
lower-level details like computation of system Jacobians can be disconnected
from the DSL and its syntax, allowing for code-reuse between modeling packages
and languages.

In this section we define the core pieces of the IR and what they mean.

### Variables

The most fundamental part of the IR is the `Variable`. In order to mirror the
intention of solving for variables and representing function-like parameters,
we treat each instance of `Variable` as a function which is called on its
arguments using the natural syntax. Rather than having additional mechanisms
for handling constant variables and parameters, we simply represent them as
constant functions.

The `Variable` is the
context-aware single variable of the IR. Its fields are described as follows:

- `name`: the name of the `Variable`. Note that this is not necessarily
  the same as the name of the Julia variable. But this symbol itself is considered
  the core identifier of the `Variable` in the sense of equality.
- `known`: the main denotation of context, storing whether or not the value of
  the variable is known.

For example, the following code defines an independent variable `t`, a parameter
`α`, a function parameter `σ`, a variable `x` which depends on `t`, a variable
`y` with no dependents, a variable `z` which depends on `t`, `α`, and `x(t)`
and a parameters `β₁` and `β₂`.


```julia
t = Variable(:t; known = true)()  # independent variables are treated as known
α = Variable(:α; known = true)()  # parameters are known
σ = Variable(:σ; known = true)    # left uncalled, since it is used as a function
w = Variable(:w; known = false)   # unknown, left uncalled
x = Variable(:x; known = false)(t)  # unknown, depends on `t`
y = Variable(:y; known = false)()   # unknown, no dependents
z = Variable(:z; known = false)(t, α, x)  # unknown, multiple arguments
β₁ = Variable(:β, 1; known = true)() # with index 1
β₂ = Variable(:β, 2; known = true)() # with index 2

expr = β₁ * x + y^α + σ(3) * (z - t) - β₂ * w(t - 1)
```

We can rewrite this more concisely using macros. Note the difference between
including and excluding empty parentheses. When in call format, variables are
aliased to the given call, allowing implicit use of dependents for convenience.

```julia
@parameters t α σ(..) β[1:2]
@variables w(..) x(t) y() z(t, α, x)

expr = β₁* x + y^α + σ(3) * (z - t) - β₂ * w(t - 1)
```

Note that `@parameters` and `@variables` implicitly add `()` to values that
are not given a call. The former specifies the values as known, while the
latter specifies it as unknown. `(..)` signifies that the value should be
left uncalled.

### Constants

`Constant` is a simple wrapper type to store numerical Julia constants.

### Operations

Operations are the basic composition of variables and puts together the pieces
with a function.

### Equations

Equations are stored using the `Equation` datatype. Given expressions for the
left-hand and right-hand sides, an equation is constructed as `Equation(lhs, rhs)`,
or equivalently `lhs ~ rhs`.

### Differentials

A `Differential` denotes the derivative with respect to a given variable. It can
be expanded via `expand_derivatives`, which symbolically differentiates
expressions recursively and cancels out appropriate constant variables.

### Systems

A system is a collection of operations with expanded context. While different
systems can have different constructors and interpretations, the general
structure is as follows:

- `eqs` is the first argument which is an array of `Operation` which describe
  the system of equations.
- Name to subtype mappings: these describe how variable `subtype`s are mapped
  to the contexts of the system. For example, for a differential equation,
  the variable corresponds to given subtypes and then the `eqs` can
  be analyzed knowing what the state variables are.
- Variable names which do not fall into one of the system's core subtypes are
  treated as intermediates which can be used for holding subcalculations and
  other pieces like that.

### Transformations

Transformation functions send IR objects to like IR objects. These utilize the
contextual information in a given `Operation`/`System` to build another
`Operation`/`System`.

## Details

### Function Registration

A function is registered into the operation system via:

```julia
@register f(x)
@register g(x,y)
```

etc. where each macro call registers the function with the given signature. This
will cause operations to stop recursing at this function, building `Operation(g,args)`
nodes into the graph instead of tracing calls of `g` itself into `Operation`s.

### Adding Derivatives

There is a large amount of derivatives pre-defined by
[DiffRules.jl](https://github.com/JuliaDiff/DiffRules.jl). Note that `Expression`
types are defined as `<:Real`, and thus any functions which allow the use of real
numbers can automatically be traced by the derivative mechanism. Thus for example:

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

where `i` means that it's the derivative of the `i`th argument. `args` is the
array of arguments, so for example if your function is `f(x,t)` then `args = [x,t]`.
You should return an `Operation` for the derivative of your function.

For example, `sin(t)`'s derivative (by `t`) is given by the following:

```julia
ModelingToolkit.derivative(::typeof(sin), args::NTuple{1,Any}, ::Val{1}) = cos(args[1])
```

### Macro-free Usage

Given the insistence on being programming friendly, all of the functionality
is accessible via a function-based interface. This means that all macros are
syntactic sugar in some form. For example, the variable construction:

```julia
@parameters t σ ρ β[1:3]
@variables x(t) y(t) z(t)
@derivatives D'~t
```

is syntactic sugar for:

```julia
t = Variable(:t; known = true)()
σ = Variable(:σ; known = true)
ρ = Variable(:ρ; known = true)()
β = [Variable(:β, i; known = true)() for i in 1:3]
x = Variable(:x)(t)
y = Variable(:y)(t)
z = Variable(:z)(t)
D = Differential(t)
```

### Intermediate Calculations

The system building functions can handle intermediate calculations. For example,

```julia
@variables x y z
@parameters σ ρ β
a = y - x
eqs = [0 ~ σ*a,
       0 ~ x*(ρ-z)-y,
       0 ~ x*y - β*z]
ns = NonlinearSystem(eqs, [x,y,z])
nlsys_func = generate_function(ns, [x,y,z], [σ,ρ,β])
```

expands to:

```julia
:((##365, u, p)->begin
          let (x, y, z, σ, ρ, β) = (u[1], u[2], u[3], p[1], p[2], p[3])
              ##365[1] = σ * (y - x)
              ##365[2] = x * (ρ - z) - y
              ##365[3] = x * y - β * z
          end
      end)
```

In addition, the Jacobian calculations take into account intermediate variables
to appropriately handle them.
