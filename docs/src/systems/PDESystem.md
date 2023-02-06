# PDESystem

`PDESystem` is the common symbolic PDE specification for the SciML ecosystem.
It is currently being built as a component of the ModelingToolkit ecosystem,

## Vision

The vision for the common PDE interface is that a user should only have to specify
their PDE once, mathematically, and have instant access to everything as simple
as a finite difference method with constant grid spacing, to something as complex
as a distributed multi-GPU discontinuous Galerkin method.

The key to the common PDE interface is a separation of the symbolic handling from
the numerical world. All the discretizers should not “solve” the PDE, but
instead be a conversion of the mathematical specification to a numerical problem.
Preferably, the transformation should be to another ModelingToolkit.jl `AbstractSystem`,
but in some cases this cannot be done or will not be performant, so a `SciMLProblem` is
the other choice.

These elementary problems, such as solving linear systems `Ax=b`, solving nonlinear
systems `f(x)=0`, ODEs, etc. are all defined by SciMLBase.jl, which then numerical
solvers can all target these common forms. Thus, someone who works on linear solvers
doesn't necessarily need to be working on a discontinuous Galerkin or finite element
library, but instead "linear solvers that are good for matrices A with
properties ..." which are then accessible by every other discretization method
in the common PDE interface.

Similar to the rest of the `AbstractSystem` types, transformation, and analysis
functions will allow for simplifying the PDE before solving it, and constructing
block symbolic functions like Jacobians.

## Constructors

```@docs
PDESystem
```

### Domains (WIP)

Domains are specifying by saying `indepvar in domain`, where `indepvar` is a
single or a collection of independent variables, and `domain` is the chosen
domain type. A 2-tuple can be used to indicate an `Interval`.
Thus forms for the `indepvar` can be like:

```julia
t ∈ (0.0, 1.0)
(t, x) ∈ UnitDisk()
[v, w, x, y, z] ∈ VectorUnitBall(5)
```

#### Domain Types (WIP)

  - `Interval(a,b)`: Defines the domain of an interval from `a` to `b` (requires explicit
    import from `DomainSets.jl`, but a 2-tuple can be used instead)

## `discretize` and `symbolic_discretize`

The only functions which act on a PDESystem are the following:

  - `discretize(sys,discretizer)`: produces the outputted `AbstractSystem` or
    `SciMLProblem`.
  - `symbolic_discretize(sys,discretizer)`: produces a debugging symbolic description
    of the discretized problem.

## Boundary Conditions (WIP)

## Transformations

## Analyses

## Discretizer Ecosystem

### NeuralPDE.jl: PhysicsInformedNN

[NeuralPDE.jl](https://docs.sciml.ai/NeuralPDE/stable/) defines the `PhysicsInformedNN`
discretizer which uses a [DiffEqFlux.jl](https://docs.sciml.ai/DiffEqFlux/stable/)
neural network to solve the differential equation.

### MethodOfLines.jl: MOLFiniteDifference

[MethodOfLines.jl](https://docs.sciml.ai/MethodOfLines/stable/) defines the
`MOLFiniteDifference` discretizer which performs a finite difference discretization.
Includes support for higher approximation order stencils and nonuniform grids.
