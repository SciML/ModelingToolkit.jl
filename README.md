# ModelingToolkit.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/ModelingToolkit/stable/)

[![codecov](https://codecov.io/gh/SciML/ModelingToolkit.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SciML/ModelingToolkit.jl)
[![Coverage Status](https://coveralls.io/repos/github/SciML/ModelingToolkit.jl/badge.svg?branch=master)](https://coveralls.io/github/SciML/ModelingToolkit.jl?branch=master)
[![Build Status](https://github.com/SciML/ModelingToolkit.jl/workflows/CI/badge.svg)](https://github.com/SciML/ModelingToolkit.jl/actions?query=workflow%3ACI)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

ModelingToolkit.jl is a modeling framework for high-performance symbolic-numeric computation
in scientific computing  and scientific machine learning.
It allows for users to give a high-level description of a model for
symbolic preprocessing to analyze and enhance the model. ModelingToolkit can
automatically generate fast functions for model components like Jacobians
and Hessians, along with automatically sparsifying and parallelizing the
computations. Automatic transformations, such as index reduction, can be applied
to the model to make it easier for numerical solvers to handle.

For information on using the package,
[see the stable documentation](https://docs.sciml.ai/ModelingToolkit/stable/). Use the
[in-development documentation](https://docs.sciml.ai/ModelingToolkit/dev/) for the version of
the documentation which contains the unreleased features.

## Standard Library

For a standard library of ModelingToolkit components and blocks, check out the
[ModelingToolkitStandardLibrary](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/)

## High-Level Examples

First, let's define a second order riff on the Lorenz equations, symbolically
lower it to a first order system, symbolically generate the Jacobian function
for the numerical integrator, and solve it.

```julia
using OrdinaryDiffEqDefault, ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters σ ρ β
@variables x(t) y(t) z(t)

eqs = [D(D(x)) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z]

@mtkbuild sys = ODESystem(eqs, t)

u0 = [D(x) => 2.0,
    x => 1.0,
    y => 0.0,
    z => 0.0]

p = [σ => 28.0,
    ρ => 10.0,
    β => 8 / 3]

tspan = (0.0, 100.0)
prob = ODEProblem(sys, u0, tspan, p, jac = true)
sol = solve(prob)
using Plots
plot(sol, idxs = (x, y))
```

![Lorenz2](https://user-images.githubusercontent.com/1814174/79118645-744eb580-7d5c-11ea-9c37-13c4efd585ca.png)

This automatically will have generated fast Jacobian functions, making
it more optimized than directly building a function. In addition, we can then
use ModelingToolkit to compose multiple ODE subsystems. Now, let's define two
interacting Lorenz equations and simulate the resulting Differential-Algebraic
Equation (DAE):

```julia
using DifferentialEquations, ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters σ ρ β
@variables x(t) y(t) z(t)

eqs = [D(x) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z]

@named lorenz1 = ODESystem(eqs, t)
@named lorenz2 = ODESystem(eqs, t)

@variables a(t)
@parameters γ
connections = [0 ~ lorenz1.x + lorenz2.y + a * γ]
@mtkbuild connected = ODESystem(connections, t, systems = [lorenz1, lorenz2])

u0 = [lorenz1.x => 1.0,
    lorenz1.y => 0.0,
    lorenz1.z => 0.0,
    lorenz2.x => 0.0,
    lorenz2.y => 1.0,
    lorenz2.z => 0.0,
    a => 2.0]

p = [lorenz1.σ => 10.0,
    lorenz1.ρ => 28.0,
    lorenz1.β => 8 / 3,
    lorenz2.σ => 10.0,
    lorenz2.ρ => 28.0,
    lorenz2.β => 8 / 3,
    γ => 2.0]

tspan = (0.0, 100.0)
prob = ODEProblem(connected, u0, tspan, p)
sol = solve(prob)

using Plots
plot(sol, idxs = (a, lorenz1.x, lorenz2.z))
```

![](https://user-images.githubusercontent.com/17304743/187790221-528046c3-dbdb-4853-b977-799596c147f3.png)

# Citation

If you use ModelingToolkit.jl in your research, please cite [this paper](https://arxiv.org/abs/2103.05244):

```
@misc{ma2021modelingtoolkit,
      title={ModelingToolkit: A Composable Graph Transformation System For Equation-Based Modeling},
      author={Yingbo Ma and Shashi Gowda and Ranjan Anantharaman and Chris Laughman and Viral Shah and Chris Rackauckas},
      year={2021},
      eprint={2103.05244},
      archivePrefix={arXiv},
      primaryClass={cs.MS}
}
```
