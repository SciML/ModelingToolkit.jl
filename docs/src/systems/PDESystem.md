# PDESystem

#### Note: PDESystem is still experimental and the solver ecosystem is being developed

`PDESystem` is the common symbolic PDE specification for the Julia ecosystem.
It is currently being built as a component of the ModelingToolkit ecosystem,
but it will soon be siphoned off to a PDE.jl which defines and documents the
whole common PDE interface ecosystem. For now, this portion documents the `PDESystem`
which is the core symbolic component of this interface.

## Vision

The vision for the common PDE interface is that a user should only have to specify
their PDE once, mathematically, and have instant access to everything as simple
as a finite difference method with constant grid spacing, to something as complex
as a distributed multi-GPU discrete Galerkin method.

The key to the common PDE interface is a separation of the symbolic handling from
the numerical world. All of the discretizers should not "solve" the PDE, but
instead be a conversion of the mathematical specification to a numerical problem.
These elementary problems, such as solving linear systems `Ax=b`, solving nonlinear
systems `f(x)=0`, ODEs, etc. are all defined by DiffEqBase.jl, which then numerical
solvers can all target these common forms. Thus someone who works on linear solvers
doesn't necessarily need to be working on a DG or finite element library, but
instead "linear solvers that are good for matrices A with properties ..." which
are then accessible by every other discretization method in the common PDE interface.
