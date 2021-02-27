# ModelingToolkit.jl

ModelingToolkit.jl is a modeling language for high-performance
symbolic-numeric computation in scientific computing and scientific machine learning.
It allows for users to give a high-level description of a model for
symbolic preprocessing to analyze and enhance the model. ModelingToolkit can
automatically generate fast functions for model components like Jacobians
and Hessians, along with automatically sparsifying and parallelizing the
computations. Automatic transformations, such as index reduction, can be applied
to the model to make it easier for numerical solvers to handle.

## Installation

To install ModelingToolkit.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("ModelingToolkit")
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
- Pervasive parallelism in symbolic computations and generated functions
- Core features like alias elimination and tearing of nonlinear systems for
  efficiently numerically handling large-scale systems of equations
- The ability to use the entire Symbolics.jl Computer Algebra System (CAS) as
  part of the modeling process.
- Extendability: the whole system is written in pure Julia, so adding new
  functions, simplification rules, and model transformations has no barrier.

For information on how to use the Symbolics.jl CAS system that ModelingToolkit.jl
is built on, consult the [Symbolics.jl documentation](https://github.com/JuliaSymbolics/Symbolics.jl)

### Equation Types

- Ordinary differential equations
- Stochastic differential equations
- Partial differential equations
- Nonlinear systems
- Optimization problems
- Optimal Control

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
- [CellMLToolkit.jl](https://github.com/SciML/CellMLToolkit.jl)
- [SbmlInterface.jl](https://github.com/paulflang/SbmlInterface.jl)
- [ReactionMechanismSimulator.jl]()
- [ReactionNetworkImporters.jl]()
- [NeuralPDE.jl]()
- [StructuralTransformations.jl]()

## Compatible Numerical Solvers



- [GalacticOptim.jl](https://github.com/SciML/GalacticOptim.jl)
