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

Let's define the Lorenz equations for numerically solving with DifferentialEquations.jl,
but tell the symbolic system to automatically generate code for efficiently
handling the sparse Jacobian.

```julia
using ModelingToolkit

@parameters t σ ρ β
@variables x(t) y(t) z(t)
@derivatives D'~t

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]

sys = ODESystem(eqs)

u0 = [x => 1.0,
      y => 0.0,
      z => 0.0]

p  = [σ => 28.0,
      ρ => 10.0,
      β => 8/3]

tspan = (0.0,100.0)      
prob = ODEProblem(sys,u0,tspan,p,jac=true,sparse=true)
```

This automatically will have generated fast sparse Jacobian functions, making
it more optimized than directly building a function. In addition, we can then
use ModelingToolkit to compose multiple ODE subsystems. Let's define two
interacting Lorenz equations:

```julia
lorenz1 = ODESystem(eqs,name=:lorenz1)
lorenz2 = ODESystem(eqs,name=:lorenz2)

@variables α
@parameters γ
connections = [0 ~ lorenz1.x + lorenz2.y + sin(α*γ)]
connected = ODESystem(connections,[α],[γ],systems=[lorenz1,lorenz2])
```

which is now a differential-algebraic equation (DAE) of 7 variables which has
two independent Lorenz systems and an algebraic equation that determines `α`
such that an implicit constraint holds. We can then define the resulting
`ODEProblem` and send it over to DifferentialEquations.jl.
