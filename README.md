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

## High Level Examples

First let's define a second order riff on the Lorenz equations, symbolically
lower it to a first order system, symbolically generate the Jacobian function
for the numerical integrator, and solve it.

```julia
using ModelingToolkit, OrdinaryDiffEq

@parameters t σ ρ β
@variables x(t) y(t) z(t)
@derivatives D'~t

eqs = [D(D(x)) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]

sys = ODESystem(eqs)
sys = ode_order_lowering(sys)

u0 = [D(x) => 2.0,
      x => 1.0,
      y => 0.0,
      z => 0.0]

p  = [σ => 28.0,
      ρ => 10.0,
      β => 8/3]

tspan = (0.0,100.0)
prob = ODEProblem(sys,u0,tspan,p,jac=true)
sol = solve(prob,Tsit5())
using Plots; plot(sol,vars=(:x,:y))
```

![Lorenz2](https://user-images.githubusercontent.com/1814174/79118645-744eb580-7d5c-11ea-9c37-13c4efd585ca.png)

This automatically will have generated fast Jacobian functions, making
it more optimized than directly building a function. In addition, we can then
use ModelingToolkit to compose multiple ODE subsystems. Now let's define two
interacting Lorenz equations and simulate the resulting Differential-Algebriac
Equation (DAE):

```julia
@parameters t σ ρ β
@variables x(t) y(t) z(t)
@derivatives D'~t

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]

lorenz1 = ODESystem(eqs,name=:lorenz1)
lorenz2 = ODESystem(eqs,name=:lorenz2)

@variables α
@parameters γ
connections = [0 ~ lorenz1.x + lorenz2.y + α*γ]
connected = ODESystem(connections,t,[α],[γ],systems=[lorenz1,lorenz2])

u0 = [lorenz1.x => 1.0,
      lorenz1.y => 0.0,
      lorenz1.z => 0.0,
	  lorenz2.x => 0.0,
	  lorenz2.y => 1.0,
	  lorenz2.z => 0.0,
	  α => 2.0]

p  = [lorenz1.σ => 10.0,
      lorenz1.ρ => 28.0,
      lorenz1.β => 8/3,
	  lorenz2.σ => 10.0,
	  lorenz2.ρ => 28.0,
	  lorenz2.β => 8/3,
	  γ => 2.0]

tspan = (0.0,100.0)
prob = ODEProblem(connected,u0,tspan,p)
sol = solve(prob,Rodas5())

using Plots; plot(sol,vars=(:α,Symbol(lorenz1.x),Symbol(lorenz2.y)))
```

![](https://user-images.githubusercontent.com/1814174/79122361-6fdaca80-7d65-11ea-87fd-0f6c4a85cd0d.png)
