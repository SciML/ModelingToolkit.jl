# ModelingToolkit.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/ModelingToolkit/stable/)

[![codecov](https://codecov.io/gh/SciML/ModelingToolkit.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SciML/ModelingToolkit.jl)
[![Coverage Status](https://coveralls.io/repos/github/SciML/ModelingToolkit.jl/badge.svg?branch=master)](https://coveralls.io/github/SciML/ModelingToolkit.jl?branch=master)
[![Build Status](https://github.com/SciML/ModelingToolkit.jl/workflows/CI/badge.svg)](https://github.com/SciML/ModelingToolkit.jl/actions?query=workflow%3ACI)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

ModelingToolkit.jl is a modeling framework for high-performance symbolic-numeric computation
in scientific computing and scientific machine learning.
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
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

# Defines a ModelingToolkit `System` model.
@parameters σ ρ β
@variables x(t) y(t) z(t)
eqs = [
    D(D(x)) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z
]
@mtkcompile sys = System(eqs, t)

# Simulate the model for a specific condition (initial condition and parameter values).
using OrdinaryDiffEqDefault
sim_cond = [
    D(x) => 2.0,
    x => 1.0,
    y => 0.0,
    z => 0.0,
    σ => 28.0,
    ρ => 10.0,
    β => 8 / 3
]
tend = 100.0
prob = ODEProblem(sys, sim_cond, tend; jac = true)
sol = solve(prob)

# Plot the solution in phase-space.
using Plots
plot(sol, idxs = (x, y))
```

![Lorenz2](https://github.com/user-attachments/assets/e82fb2ce-97b7-4f56-b272-85653c88bdb3)

This will have automatically generated fast Jacobian functions, making
it more optimized than directly building a function. In addition, we can then
use ModelingToolkit to compose multiple ODE subsystems. Now, let's define two
interacting Lorenz equations and simulate the resulting Differential-Algebraic
Equation (DAE):

```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

# Defines two lorenz system models.
eqs = [
    D(x) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z
]
@named lorenz1 = System(eqs, t)
@named lorenz2 = System(eqs, t)

# Connect the two models, creating a single model.
@variables a(t)
@parameters γ
connections = [0 ~ lorenz1.x + lorenz2.y + a * γ]
@mtkcompile connected_lorenz = System(connections, t; systems = [lorenz1, lorenz2])

# Simulate the model for a specific condition (initial condition and parameter values).
using OrdinaryDiffEqDefault
sim_cond = [
    lorenz1.x => 1.0,
    lorenz1.y => 0.0,
    lorenz1.z => 0.0,
    lorenz2.x => 0.0,
    lorenz2.z => 0.0,
    a => 2.0,
    lorenz1.σ => 10.0,
    lorenz1.ρ => 28.0,
    lorenz1.β => 8 / 3,
    lorenz2.σ => 10.0,
    lorenz2.ρ => 28.0,
    lorenz2.β => 8 / 3,
    γ => 2.0
]
tend = 100.0
prob = ODEProblem(connected_lorenz, sim_cond, tend)
sol = solve(prob)

# Plot the solution in phase-space.
using Plots
plot(sol, idxs = (a, lorenz1.x, lorenz2.z))
```

![LorenzConnected](https://github.com/user-attachments/assets/ef65d812-c10e-42c6-945d-e61515e3b6a1)

## Citation

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
