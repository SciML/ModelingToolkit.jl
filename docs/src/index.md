# ModelingToolkit.jl

ModelingToolkit.jl is an intermediate representation (IR) of computational graphs
for scientific computing problems. Its purpose is to be a common target for
modeling domain-specific languages (DSLs) in order to allow for a common
platform for model inspection and transformation. It uses a tagged variable IR
in order to allow specification of complex models and allow for transformations
of models. It has ways to plug into its function registration and derivative
system so that way it can interact nicely with user-defined routines. Together,
this is an abstract form of a scientific model that is easy for humans to
generate but also easy for programs to manipulate.

## Package Overview

ModelingToolkit has 4 layers:

1. The high level. This level has syntactic sugar for easily generating
   ModelingToolkit models. It can be used directly like a DSL for advanced
   users who want flexibility, and for easily generating DSLs.
2. The tracing level. This level allows for easily generating ModelingToolkit
   models directly from Julia code by tracing existing functions and turning
   them into ModelingToolkit IR and `AbstractSystem`s for symbolic manipulation.
3. The `AbstractSystem` level. This is the level where content-dependent functionality
   is added, where models such an ordinary differential equation are represented.
   At the system level there are *transformations* which take one system to
   another, and *targets* which output code for numerical solvers.
4. The IR level, also referred to as the direct level. At this level, one
   directly acts on arrays of `Equation`, `Operation` and `Variable` types to
   generate functions.

Each level above is built on the level below, giving more context to allow for
more automation. For example, the system level allows for automatically generating
fast multithreaded sparse Jacobian functions of an `ODESystem`, which is just
calling the sparsity functions and the multithreading capabilities of
`build_function` at the IR level.
