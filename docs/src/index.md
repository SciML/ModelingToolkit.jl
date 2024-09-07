# ModelingToolkit.jl: High-Performance Symbolic-Numeric Equation-Based Modeling

ModelingToolkit.jl is a modeling language for high-performance
symbolic-numeric computation in scientific computing and scientific machine learning.
It then mixes ideas from symbolic computational algebra systems with
causal and acausal equation-based modeling frameworks to give an extendable and
parallel modeling system. It allows for users to give a high-level description of
a model for symbolic preprocessing to analyze and enhance the model. Automatic
symbolic transformations, such as index reduction of differential-algebraic equations,
make it possible to solve equations that are impossible to solve
with a purely numeric-based technique.

## Installation

To install ModelingToolkit.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("ModelingToolkit")
```

## Citation

If you use ModelingToolkit in your work, please cite the following:

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

## Feature Summary

ModelingToolkit.jl is a symbolic-numeric modeling package. Thus it combines some
of the features from symbolic computing packages like SymPy or Mathematica with
the ideas of equation-based modeling systems like the causal Simulink and the
acausal Modelica. It bridges the gap between many different kinds of equations,
allowing one to quickly and easily transform systems of DAEs into optimization
problems, or vice-versa, and then simplify and parallelize the resulting expressions
before generating code.

### Feature List

  - Causal and acausal modeling (Simulink/Modelica)
  - Automated model transformation, simplification, and composition
  - Automatic conversion of numerical models into symbolic models
  - Composition of models through the components, a lazy connection system, and
    tools for expanding/flattening
  - Pervasive parallelism in symbolic computations and generated functions
  - Transformations like alias elimination and tearing of nonlinear systems for
    efficiently numerically handling large-scale systems of equations
  - The ability to use the entire Symbolics.jl Computer Algebra System (CAS) as
    part of the modeling process.
  - Import models from common formats like SBML, CellML, BioNetGen, and more.
  - Extensibility: the whole system is written in pure Julia, so adding new
    functions, simplification rules, and model transformations has no barrier.

For information on how to use the Symbolics.jl CAS system that ModelingToolkit.jl
is built on, consult the
[Symbolics.jl documentation](https://docs.sciml.ai/Symbolics/stable/)

### Equation Types

  - Ordinary differential equations
  - Stochastic differential equations
  - Partial differential equations
  - Nonlinear systems
  - Optimization problems
  - Continuous-Time Markov Chains
  - Chemical Reactions (via [Catalyst.jl](https://docs.sciml.ai/Catalyst/stable/))
  - Nonlinear Optimal Control

## Standard Library

For quick development, ModelingToolkit.jl includes
[ModelingToolkitStandardLibrary.jl](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/),
a standard library of prebuilt components for the ModelingToolkit ecosystem.

## Model Import Formats

  - [CellMLToolkit.jl](https://docs.sciml.ai/CellMLToolkit/stable/): Import [CellML](https://www.cellml.org/) models into ModelingToolkit
    
      + Repository of more than a thousand pre-made models
      + Focus on biomedical models in areas such as: Calcium Dynamics,
        Cardiovascular Circulation, Cell Cycle, Cell Migration, Circadian Rhythms,
        Electrophysiology, Endocrine, Excitation-Contraction Coupling, Gene Regulation,
        Hepatology, Immunology, Ion Transport, Mechanical Constitutive Laws,
        Metabolism, Myofilament Mechanics, Neurobiology, pH Regulation, PKPD,
        Protein Modules, Signal Transduction, and Synthetic Biology.

  - [SBMLToolkit.jl](https://docs.sciml.ai/SBMLToolkit/stable/): Import [SBML](http://sbml.org/) models into ModelingToolkit
    
      + Uses the robust libsbml library for parsing and transforming the SBML
  - [ReactionNetworkImporters.jl](https://docs.sciml.ai/ReactionNetworkImporters/stable/): Import various models into ModelingToolkit
    
      + Supports the BioNetGen `.net` file
      + Supports importing networks specified by stoichiometric matrices

## Extension Libraries

Because ModelingToolkit.jl is the core foundation of an equation-based modeling
ecosystem, there is a large set of libraries adding features to this system.
Below is an incomplete list of extension libraries one may want to be aware of:

  - [Catalyst.jl](https://docs.sciml.ai/Catalyst/stable/): Symbolic representations
    of chemical reactions
    
      + Symbolically build and represent large systems of chemical reactions
      + Generate code for ODEs, SDEs, continuous-time Markov Chains, and more
      + Simulate the models using the SciML ecosystem with O(1) Gillespie methods

  - [DataDrivenDiffEq.jl](https://docs.sciml.ai/DataDrivenDiffEq/stable/): Automatic
    identification of equations from data
    
      + Automated construction of ODEs and DAEs from data
      + Representations of Koopman operators and Dynamic Mode Decomposition (DMD)
  - [MomentClosure.jl](https://augustinas1.github.io/MomentClosure.jl/dev/): Automatic
    transformation of ReactionSystems into deterministic systems
    
      + Generates ODESystems for the moment closures
      + Allows for geometrically-distributed random reaction rates
  - [ReactionMechanismSimulator.jl](https://github.com/ReactionMechanismGenerator/ReactionMechanismSimulator.jl):
    Simulating and analyzing large chemical reaction mechanisms
    
      + Ideal gas and dilute liquid phases.
      + Constant T and P and constant V adiabatic ideal gas reactors.
      + Constant T and V dilute liquid reactors.
      + Diffusion limited rates. Sensitivity analysis for all reactors.
      + Flux diagrams with molecular images (if molecular information is provided).
  - [NumCME.jl](https://github.com/voduchuy/NumCME.jl): High-performance simulation of chemical master equations (CME)
    
      + Transient solution of the CME
      + Dynamic state spaces
      + Accepts reaction systems defined using Catalyst.jl DSL.
  - [FiniteStateProjection.jl](https://github.com/SciML/FiniteStateProjection.jl): High-performance simulation of
    chemical master equations (CME) via finite state projections
    
      + Accepts reaction systems defined using Catalyst.jl DSL.

## Compatible Numerical Solvers

All of the symbolic systems have a direct conversion to a numerical system, which
can then be handled through the SciML interfaces. For example, after building a
model and performing symbolic manipulations, an `ODESystem` can be converted into
an `ODEProblem` to then be solved by a numerical ODE solver. Below is a list of
the solver libraries which are the numerical targets of the ModelingToolkit
system:

  - [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/)
    
      + Multi-package interface of high performance numerical solvers for `ODESystem`,
        `SDESystem`, and `JumpSystem`

  - [NonlinearSolve.jl](https://docs.sciml.ai/NonlinearSolve/stable/)
    
      + High performance numerical solving of `NonlinearSystem`
  - [Optimization.jl](https://docs.sciml.ai/Optimization/stable/)
    
      + Multi-package interface for numerical solving `OptimizationSystem`
  - [NeuralPDE.jl](https://docs.sciml.ai/NeuralPDE/stable/)
    
      + Physics-Informed Neural Network (PINN) training on `PDESystem`
  - [MethodOfLines.jl](https://docs.sciml.ai/MethodOfLines/stable/)
    
      + Automated finite difference method (FDM) discretization of `PDESystem`

## Contributing

  - Please refer to the
    [SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://github.com/SciML/ColPrac/blob/master/README.md)
    for guidance on PRs, issues, and other matters relating to contributing to SciML.

  - See the [SciML Style Guide](https://github.com/SciML/SciMLStyle) for common coding practices and other style decisions.
  - There are a few community forums:
    
      + The #diffeq-bridged and #sciml-bridged channels in the
        [Julia Slack](https://julialang.org/slack/)
      + The #diffeq-bridged and #sciml-bridged channels in the
        [Julia Zulip](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
      + On the [Julia Discourse forums](https://discourse.julialang.org)
      + See also [SciML Community page](https://sciml.ai/community/)

## Reproducibility

```@raw html
<details><summary>The documentation of this SciML package was built using these direct dependencies,</summary>
```

```@example
using Pkg # hide
Pkg.status() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>and using this machine and Julia version.</summary>
```

```@example
using InteractiveUtils # hide
versioninfo() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>
```

```@example
using Pkg # hide
Pkg.status(; mode = PKGMODE_MANIFEST) # hide
```

```@raw html
</details>
```

```@eval
using TOML
using Markdown
version = TOML.parse(read("../../Project.toml", String))["version"]
name = TOML.parse(read("../../Project.toml", String))["name"]
link_manifest = "https://github.com/SciML/" * name * ".jl/tree/gh-pages/v" * version *
                "/assets/Manifest.toml"
link_project = "https://github.com/SciML/" * name * ".jl/tree/gh-pages/v" * version *
               "/assets/Project.toml"
Markdown.parse("""You can also download the
[manifest]($link_manifest)
file and the
[project]($link_project)
file.
""")
```
