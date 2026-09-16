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

Dependent variables may use the standard ModelingToolkit input and output metadata.
The declarations are available through `inputs(sys)` and `outputs(sys)` in dependent-variable
declaration order. If a variable is marked as both an input and an output, it is reported as an
input. These roles describe the symbolic PDE interface; support for discretizing them depends on
the selected PDE discretizer.

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

## Solution Interface

Whatever the discretizer, `solve(prob, alg)` on the problem returned by `discretize`
gives back a solution expressed in the `PDESystem`'s own variables: a
`PDETimeSeriesSolution` when the system has a time variable and a `PDENoTimeSolution`
otherwise (both from SciMLBase). Every discretizer indexes and evaluates them the same
way:

  - `sol[u(t, x)]` is the dependent variable `u` on the discretization grid (or, for a
    mesh-free method, on its evaluation grid), as an array with one axis per argument of
    `u`.
  - `sol[x]` is the grid of the independent variable `x`; for a time-dependent solution
    `sol[t]` (also `sol.t`) holds the saved times.
  - `sol(t, x; dv = u(t, x))` evaluates `u` at arbitrary points: numbers or ranges, one
    per independent variable, interpolated on a grid-based discretization and evaluated
    directly by a mesh-free one. Without `dv` it returns the values of every dependent
    variable.
  - `sol.original_sol` is the underlying `ODESolution`, `OptimizationSolution`, or other
    solution of the discretized problem, for anything the wrapper does not expose.

[MethodOfLines.jl](https://docs.sciml.ai/MethodOfLines/stable/solutions/) and
[NeuralPDE.jl](https://docs.sciml.ai/NeuralPDE/stable/) implement this interface; see the
[PDEBase.jl developer documentation](https://docs.sciml.ai/PDEBase/stable/interface/) for
what a new discretizer has to define.

## Boundary Conditions (WIP)

## Transformations

## Analyses

## Discretizer Ecosystem

### NeuralPDE.jl: PhysicsInformedNN

[NeuralPDE.jl](https://docs.sciml.ai/NeuralPDE/stable/) defines the `PhysicsInformedNN`
discretizer, a physics-informed neural network: the dependent variables are represented by
[Lux.jl](https://lux.csail.mit.edu/) networks whose parameters become the unknowns of an
`OptimizationProblem` minimizing the residuals on collocation points.

### MethodOfLines.jl: MOLFiniteDifference

[MethodOfLines.jl](https://docs.sciml.ai/MethodOfLines/stable/) defines the
`MOLFiniteDifference` discretizer which performs a finite difference discretization.
Includes support for higher approximation order stencils and nonuniform grids.
