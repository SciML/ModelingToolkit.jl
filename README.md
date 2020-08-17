# ModelingToolkit.jl

[![Build Status](https://travis-ci.com/SciML/ModelingToolkit.jl.svg?branch=master)](https://travis-ci.com/SciML/ModelingToolkit.jl)
[![Coverage Status](https://coveralls.io/repos/SciML/ModelingToolkit.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JuliaDiffEq/ModelingToolkit.jl?branch=master)
[![codecov.io](http://codecov.io/github/SciML/ModelingToolkit.jl/coverage.svg?branch=master)](http://codecov.io/github/SciML/ModelingToolkit.jl?branch=master)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](http://mtk.sciml.ai/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](http://mtk.sciml.ai/dev/)

ModelingToolkit.jl is a modeling framework for high-performance symbolic-numeric computation
in scientific computing  and scientific machine learning.
It allows for users to give a high-level description of a model for
symbolic preprocessing to analyze and enhance the model. ModelingToolkit can
automatically generate fast functions for model components like Jacobians
and Hessians, along with automatically sparsifying and parallelizing the
computations. Automatic transformations, such as index reduction, can be applied
to the model to make it easier for numerical solvers to handle.

For information on using the package,
[see the stable documentation](https://mtk.sciml.ai/stable/). Use the
[in-development documentation](https://mtk.sciml.ai/dev/) for the version of
the documentation which contains the unreleased features.

## High-Level Examples

First, let's define a second order riff on the Lorenz equations, symbolically
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

p  = [σ => 10.0,
      ρ => 28.0,
      β => 8/3]

tspan = (0.0,100.0)
prob = ODEProblem(sys,u0,tspan,p,jac=true)
sol = solve(prob,Tsit5())
using Plots; plot(sol,vars=(x,y))
```

![Lorenz2](https://user-images.githubusercontent.com/1814174/79118645-744eb580-7d5c-11ea-9c37-13c4efd585ca.png)

This automatically will have generated fast Jacobian functions, making
it more optimized than directly building a function. In addition, we can then
use ModelingToolkit to compose multiple ODE subsystems. Now, let's define two
interacting Lorenz equations and simulate the resulting Differential-Algebraic
Equation (DAE):

```julia
using ModelingToolkit, OrdinaryDiffEq

@parameters t σ ρ β
@variables x(t) y(t) z(t)
@derivatives D'~t

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]

lorenz1 = ODESystem(eqs,name=:lorenz1)
lorenz2 = ODESystem(eqs,name=:lorenz2)

@variables a
@parameters γ
connections = [0 ~ lorenz1.x + lorenz2.y + a*γ]
connected = ODESystem(connections,t,[a],[γ],systems=[lorenz1,lorenz2])

u0 = [lorenz1.x => 1.0,
      lorenz1.y => 0.0,
      lorenz1.z => 0.0,
      lorenz2.x => 0.0,
      lorenz2.y => 1.0,
      lorenz2.z => 0.0,
      a => 2.0]

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

using Plots; plot(sol,vars=(a,lorenz1.x,lorenz2.z))
```

![](https://user-images.githubusercontent.com/1814174/79229194-9e71a780-7e30-11ea-9f93-bfa762eb8dfb.png)
