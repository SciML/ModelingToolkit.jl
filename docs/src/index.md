# ModelingToolkit.jl: High-Performance Symbolic-Numeric Equation-Based Modeling

ModelingToolkit.jl is a modeling language for high-performance
symbolic-numeric computation in scientific computing and scientific machine learning.
It then mixes ideas from symbolic computational algebra systems with
causal and acausal equation-based modeling frameworks to give an extendable and
parallel modeling system. It allows for users to give a high-level description of
a model for symbolic preprocessing to analyze and enhance the model. Automatic
transformations, such as index reduction, can be applied to the model before
solving in order to make it easily handle equations would could not be solved
when modeled without symbolic intervention.

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
- Extendability: the whole system is written in pure Julia, so adding new
  functions, simplification rules, and model transformations has no barrier.

For information on how to use the Symbolics.jl CAS system that ModelingToolkit.jl
is built on, consult the
[Symbolics.jl documentation](https://github.com/JuliaSymbolics/Symbolics.jl)

### Equation Types

- Ordinary differential equations
- Stochastic differential equations
- Partial differential equations
- Nonlinear systems
- Optimization problems
- Continuous-Time Markov Chains
- Chemical Reactions
- Nonlinear Optimal Control

## Model Import Formats

- [CellMLToolkit.jl](https://github.com/SciML/CellMLToolkit.jl): Import [CellML](https://www.cellml.org/) models into ModelingToolkit
    - Repository of more than a thousand pre-made models
    - Focus on biomedical models in areas such as: Calcium Dynamics,
      Cardiovascular Circulation, Cell Cycle, Cell Migration, Circadian Rhythms,
      Electrophysiology, Endocrine, Excitation-Contraction Coupling, Gene Regulation,
      Hepatology, Immunology, Ion Transport, Mechanical Constitutive Laws,
      Metabolism, Myofilament Mechanics, Neurobiology, pH Regulation, PKPD,
      Protein Modules, Signal Transduction, and Synthetic Biology.
- [SbmlInterface.jl](https://github.com/paulflang/SbmlInterface.jl): Import [SBML](http://sbml.org/Main_Page) models into ModelingToolkit
    - Uses the robust libsbml library for parsing and transforming the SBML
- [ReactionNetworkImporters.jl](https://github.com/isaacsas/ReactionNetworkImporters.jl): Import various models into ModelingToolkit
    - Supports the BioNetGen `.net` file
    - Supports importing networks specified by stoichiometric matrices

## Extension Libraries

Because ModelingToolkit.jl is the core foundation of a equation-based modeling
ecosystem, there is a large set of libraries adding features to this system.
Below is an incomplete list of extension libraries one may want to be aware of:

- [Catalyst.jl](https://github.com/SciML/Catalyst.jl): Symbolic representations of chemical reactions
    - Symbolically build and represent large systems of chemical reactions
    - Generate code for ODEs, SDEs, continuous-time Markov Chains, and more
    - Simulate the models using the SciML ecosystem with O(1) Gillespie methods
- [DataDrivenDiffEq.jl](https://github.com/SciML/DataDrivenDiffEq.jl): Automatic identification of equations from data
    - Automated construction of ODEs and DAEs from data
    - Representations of Koopman operators and Dynamic Mode Decomposition (DMD)
- [MomentClosure.jl](https://github.com/augustinas1/MomentClosure.jl): Automatic transformation of ReactionSystems into deterministic systems
    - Generates ODESystems for the moment closures
    - Allows for geometrically-distributed random reaction rates
- [ReactionMechanismSimulator.jl](https://github.com/ReactionMechanismGenerator/ReactionMechanismSimulator.jl): simulating and analyzing large chemical reaction mechanisms
    - Ideal gas and dilute liquid phases.
    - Constant T and P and constant V adiabatic ideal gas reactors.
    - Constant T and V dilute liquid reactors.
    - Diffusion limited rates. Sensitivity analysis for all reactors.
    - Flux diagrams with molecular images (if molecular information is provided).

## Compatible Numerical Solvers

All of the symbolic systems have a direct conversion to a numerical system which
can then be handled through the SciML interfaces. For example, after building a
model and performing symbolic manipulations, an `ODESystem` can be converted into
an `ODEProblem` to then be solved by a numerical ODE solver. Below is a list of
the solver libraries which are the numerical targets of the ModelingToolkit
system:

- [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/)
    - Multi-package interface of high performance numerical solvers for `ODESystem`, `SDESystem`, and `JumpSystem`
- [NonlinearSolve.jl](https://github.com/JuliaComputing/NonlinearSolve.jl)
    - High performance numerical solving of `NonlinearSystem`
- [GalacticOptim.jl](https://github.com/SciML/GalacticOptim.jl)
    - Multi-package interface for numerical solving `OptimizationSystem`
- [NeuralPDE.jl](https://github.com/SciML/NeuralPDE.jl)
    - Physics-Informed Neural Network (PINN) training on `PDESystem`
- [DiffEqOperators.jl](https://github.com/SciML/DiffEqOperators.jl)
    - Automated finite difference method (FDM) discretization of `PDESystem`

## Contributing

- Please refer to the [SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://github.com/SciML/ColPrac/blob/master/README.md) for guidance on PRs, issues, and other matters relating to contributing to ModelingToolkit.
- There are a few community forums:
    - The #diffeq-bridged channel in the [Julia Slack](https://julialang.org/slack/)
    - [JuliaDiffEq](https://gitter.im/JuliaDiffEq/Lobby) on Gitter
    - On the Julia Discourse forums (look for the [modelingtoolkit tag](https://discourse.julialang.org/tag/modelingtoolkit)
    - See also [SciML Community page](https://sciml.ai/community/)
