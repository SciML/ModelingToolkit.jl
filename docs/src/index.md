# ModelingToolkit.jl

ModelingToolkit.jl is a modeling language for scientific computing problems.
It allows for users to give a high-level description of a model for
symbolic preprocessing to analyze and enhance the model. ModelingToolkit can
automatically generate fast functions for model components like Jacobians
and Hessians, along with automatically sparsifying and parallelizing the
computations. Automatic transformations, such as index reduction, can be applied
to the model to make it easier for numerical solvers to handle.

## Package Overview

ModelingToolkit has 3 layers:

1. The model definition level. This is a high level of syntactic sugar for
   easily generating ModelingToolkit models. It can be used directly like a DSL
   for advanced users who want a lot of flexibility in a modeling language.
   Additionally, automatic tracing functionality allows for easily generating
   ModelingToolkit models directly from Julia code.
2. The `AbstractSystem` level. This is the level where content-dependent functionality
   is added, where models such an ordinary differential equation are represented.
   At the system level there are *transformations* which take one system to
   another, and *targets* which output code for numerical solvers.
3. The IR level, also referred to as the direct level. At this level, one
   directly acts on arrays of `Equation`, `Operation` and `Variable` types to
   generate functions.

Each level above is built on the level below, giving more context to allow for
more automation. For example, the system level allows for automatically generating
fast multithreaded sparse Jacobian functions of an `ODESystem`, which is just
calling the sparsity functions and the multithreading capabilities of
`build_function` at the IR level.
