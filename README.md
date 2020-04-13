# ModelingToolkit.jl

[![Build Status](https://travis-ci.org/SciML/ModelingToolkit.jl.svg?branch=master)](https://travis-ci.org/SciML/ModelingToolkit.jl)
[![Coverage Status](https://coveralls.io/repos/SciML/ModelingToolkit.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JuliaDiffEq/ModelingToolkit.jl?branch=master)
[![codecov.io](http://codecov.io/github/SciML/ModelingToolkit.jl/coverage.svg?branch=master)](http://codecov.io/github/SciML/ModelingToolkit.jl?branch=master)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](http://mtk.sciml.ai/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](http://mtk.sciml.ai/dev/)

ModelingToolkit.jl is an intermediate representation (IR) of computational graphs
for scientific computing problems. Its purpose is to be a common target for
modeling DSLs in order to allow for a common platform for model inspection and
transformation. It uses a tagged variable IR in order to allow specification of
complex models and allow for transformations of models. It has ways to plug into
its function registration and derivative system so that way it can interact
nicely with user-defined routines. Together, this is an abstract form of a
scientific model that is easy for humans to generate but also easy for programs
to manipulate.

For information on using the package,
[see the stable documentation](https://mtk.sciml.ai/stable/). Use the
[in-development documentation](https://mtk.sciml.ai/dev/) for the version of
the documentation which contains the un-released features.

## High Level Example

First we define some variables. In a differential equation
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
de = ODESystem(eqs, t, [x,y,z], [σ,ρ,β])
```

where we tell it the variable types and ordering in the first version, or let it
automatically determine the variable types in the second version.
This can then generate the function. For example, we can see the
generated code via:

```julia
myode_oop = generate_function(de,linenumbers=false)[1] # first one is the out-of-place function

#=
:((u, p, t)->begin
          if u isa Array
              return @inbounds(begin
                          let (x, y, z, σ, ρ, β) = (u[1], u[2], u[3], p[1], p[2], p[3])
                              [σ * (y - x), x * (ρ - z) - y, x * y - β * z]
                          end
                      end)
          else
              X = @inbounds(begin
                          let (x, y, z, σ, ρ, β) = (u[1], u[2], u[3], p[1], p[2], p[3])
                              (σ * (y - x), x * (ρ - z) - y, x * y - β * z)
                          end
                      end)
          end
          T = promote_type(map(typeof, X)...)
          map(T, X)
          construct = if u isa ModelingToolkit.StaticArrays.StaticArray
                  ModelingToolkit.StaticArrays.similar_type(typeof(u), eltype(X))
              else
                  x->begin
                          convert(typeof(u), x)
                      end
              end
          construct(X)
      end)
=#

myode_iip = generate_function(de,linenumbers=false)[2] # second one is the in-place function

#=
:((var"##MTIIPVar#793", u, p, t)->begin
          @inbounds begin
                  @inbounds begin
                          let (x, y, z, σ, ρ, β) = (u[1], u[2], u[3], p[1], p[2], p[3])
                              var"##MTIIPVar#793"[1] = σ * (y - x)
                              var"##MTIIPVar#793"[2] = x * (ρ - z) - y
                              var"##MTIIPVar#793"[3] = x * y - β * z
                          end
                      end
              end
          nothing
      end)
=#
```

or directly get the generated ODE function via:

```julia
f = ODEFunction(de, [x,y,z], [σ,ρ,β])
```

Here already you can see some advantages of the ModelingToolkit.jl compilation system. As an
IR to target, this output can compile to multiple different forms, including ones specific
to static arrays and in-place functions. Forms which automatically parallelize the calculations
based on internal cost models are a work-in-progress as well. This means DSLs built on top of
this as a model compiler can write domain-specific languages without having to write complex
optimized Julia function compilers.

## Tutorial

For an introductory tutorial to using ModelingToolkit.jl, please checkout
[ModelingToolkit.jl, An IR and Compiler for Scientific Models](https://tutorials.juliadiffeq.org/html/ode_extras/01-ModelingToolkit.html).
