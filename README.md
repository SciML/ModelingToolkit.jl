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

#### Warning: This repository is a work-in-progress

## Introduction by Examples

### Example: ODE

Let's build an ODE. First we define some variables. In a differential equation
system, we need to differentiate between our dependent variables, independent
variables, and parameters. Therefore we label them as follows:

```julia
using ModelingToolkit

# Define some variables
@IVar t
@DVar x(t) y(t) z(t)
@Deriv D'~t
@Param σ ρ β
```

Then we build the system:

```julia
eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]
```

Each operation builds an `Operation` type, and thus `eqs` is an array of
`Operation` and `Variable`s. This holds a tree of the full system that can be
analyzed by other programs. We can turn this into a `DiffEqSystem` via:

```julia
de = DiffEqSystem(eqs,[t],[x,y,z],Variable[],[σ,ρ,β])
de = DiffEqSystem(eqs)
```

where we tell it the variable types and ordering in the first version, or let it
automatically determine the variable types in the second version.
This can then generate the function. For example, we can see the
generated code via:

```julia
ModelingToolkit.generate_ode_function(de)

## Which returns:
:((du, u, p, t)->begin
                x = u[1]
                y = u[2]
                z = u[3]
                σ = p[1]
                ρ = p[2]
                β = p[3]
                x_t = σ * (y - x)
                y_t = x * (ρ - z) - y
                z_t = x * y - β * z
                du[1] = x_t
                du[2] = y_t
                du[3] = z_t
            end
        end)
```

and get the generated function via:

```julia
f = ODEFunction(de)
```

### Example: Nonlinear System

We can also build nonlinear systems. Let's say we wanted to solve for the steady
state of the previous ODE. This is the nonlinear system defined by where the
derivatives are zero. We could use dependent variables for our nonlinear system
(for direct compatibility with the above ODE code), or we can use non-tagged
variables. Here we will show the latter. We write:

```julia
@Var x y z
@Param σ ρ β

# Define a nonlinear system
eqs = [0 ~ σ*(y-x),
       0 ~ x*(ρ-z)-y,
       0 ~ x*y - β*z]
ns = NonlinearSystem(eqs)
nlsys_func = ModelingToolkit.generate_nlsys_function(ns)
```

which generates:

```julia
(du, u, p)->begin  # C:\Users\Chris\.julia\v0.6\ModelingToolkit\src\systems.jl, line 51:
        begin  # C:\Users\Chris\.julia\v0.6\ModelingToolkit\src\utils.jl, line 2:
            y = u[1]
            x = u[2]
            z = u[3]
            σ = p[1]
            ρ = p[2]
            β = p[3]
            resid[1] = σ * (y - x)
            resid[2] = x * (ρ - z) - y
            resid[3] = x * y - β * z
        end
    end
```

We can use this to build a nonlinear function for use with NLsolve.jl:

```julia
f = @eval eval(nlsys_func)
# Make a closure over the parameters for for NLsolve.jl
f2 = (du,u) -> f(du,u,(10.0,26.0,2.33))
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

The most fundamental part of the IR is the `Variable`. The `Variable` is the
context-aware single variable of the IR. Its fields are described as follows:

- `name`: the name of the `Variable`. Note that this is not necessarily
  the same as the name of the Julia variable. But this symbol itself is considered
  the core identifier of the `Variable` in the sense of equality.
- `subtype`: the main denotation of context. Variables within systems
  are grouped according to their `subtype`.
- `diff`: the `Differential` object representing the quantity the variable is differentiated with respect to, or `nothing`
- `dependents`: the vector of variables on which the current variable
  is dependent. For example, `u(t,x)` has dependents `[t,x]`. Derivatives thus
  require this information in order to simplify down.

### Constants

`Constant` is a simple wrapper type to store numerical Julia constants.

### Operations

Operations are the basic composition of variables and puts together the pieces
with a function. The `~` function denotes equality between the arguments.

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
  the dependent variable corresponds to given subtypes and then the `eqs` can
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
ModelingToolkit.Derivative(::typeof(my_function),args,::Val{i})
```

where `i` means that it's the derivative of the `i`th argument. `args` is the
array of arguments, so for example if your function is `f(x,t)` then `args = [x,t]`.
You should return an `Operation` for the derivative of your function.

For example, `sin(t)`'s derivative (by `t`) is given by the following:

```julia
ModelingToolkit.Derivative(::typeof(sin),args,::Val{1}) = cos(args[1])
```

### Macro-free Usage

Given the insistence on being programming friendly, all of the functionality
is accessible via a function-based interface. This means that all macros are
syntactic sugar in some form. For example, the variable construction:

```julia
@IVar t
@DVar x(t) y(t) z(t)
@Deriv D'~t
@Param σ ρ β
```

is syntactic sugar for:

```julia
t = IndependentVariable(:t)
x = DependentVariable(:x,dependents = [t])
y = DependentVariable(:y,dependents = [t])
z = DependentVariable(:z,dependents = [t])
D = Differential(t) # Default of first derivative, Derivative(t,1)
σ = Parameter(:σ)
ρ = Parameter(:ρ)
β = Parameter(:β)
```

### Intermediate Calculations

The system building functions can handle intermediate calculations. For example,

```julia
@Var a x y z
@Param σ ρ β
eqs = [a ~ y-x,
       0 ~ σ*a,
       0 ~ x*(ρ-z)-y,
       0 ~ x*y - β*z]
ns = NonlinearSystem(eqs,[x,y,z],[σ,ρ,β])
nlsys_func = ModelingToolkit.generate_nlsys_function(ns)
```

expands to:

```julia
:((du, u, p)->begin  # C:\Users\Chris\.julia\v0.6\ModelingToolkit\src\systems.jl, line 85:
            begin  # C:\Users\Chris\.julia\v0.6\ModelingToolkit\src\utils.jl, line 2:
                x = u[1]
                y = u[2]
                z = u[3]
                σ = p[1]
                ρ = p[2]
                β = p[3]
                a = y - x
                resid[1] = σ * a
                resid[2] = x * (ρ - z) - y
                resid[3] = x * y - β * z
            end
        end)
```

In addition, the Jacobian calculations take into account intermediate variables
to appropriately handle them.
