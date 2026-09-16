```@meta
CollapsedDocStrings = true
```

# Building and solving numerical problems

Systems are numerically solved by building and solving the appropriate problem type.
Numerical solvers expect to receive functions taking a predefeined set of arguments
and returning specific values. This format of argument and return value depends on
the function and the problem. ModelingToolkit is capable of compiling and generating
code for a variety of such numerical problems.

## In-place and out-of-place problems

Every problem and function constructor takes an `iip` type parameter selecting an in-place
or out-of-place formulation. Problem constructors called without it defer the choice to
construction time using a sentinel type.

```@docs
ModelingToolkitBase.Both
```

## Dynamical systems

```@docs
SciMLBase.ODEFunction
SciMLBase.ODEProblem
SciMLBase.DAEFunction
SciMLBase.DAEProblem
ModelingToolkit.SemilinearODEFunction
ModelingToolkit.SemilinearODEProblem
SciMLBase.SDEFunction
SciMLBase.SDEProblem
SciMLBase.DDEFunction
SciMLBase.DDEProblem
SciMLBase.SDDEFunction
SciMLBase.SDDEProblem
JumpProcesses.JumpProblem
SciMLBase.BVProblem
SciMLBase.DiscreteFunction
SciMLBase.DiscreteProblem
SciMLBase.ImplicitDiscreteFunction
SciMLBase.ImplicitDiscreteProblem
```

## Linear and Nonlinear systems

```@docs
SciMLBase.NonlinearFunction
SciMLBase.NonlinearProblem
SciMLBase.AbstractNonlinearProblem(::System, ::Any)
SciMLBase.HomotopyProblem
SciMLBase.SCCNonlinearProblem
SciMLBase.NonlinearLeastSquaresProblem
SciMLBase.SteadyStateProblem
SciMLBase.IntervalNonlinearFunction
SciMLBase.IntervalNonlinearProblem
ModelingToolkit.HomotopyContinuationProblem
SciMLBase.HomotopyNonlinearFunction
SciMLBase.LinearProblem
```

## Problem type metadata

Systems can carry a `problem_type` that is forwarded to the constructed problem. This is
how discretization packages attach their own metadata to a generated problem and recover
it downstream, for example to wrap a solution in a richer solution type.

```@docs
ModelingToolkit.ProblemTypeCtx
```

## Optimization and optimal control

```@docs
SciMLBase.OptimizationFunction
SciMLBase.MultiObjectiveOptimizationFunction
SciMLBase.OptimizationProblem
SciMLBase.ODEInputFunction
ModelingToolkit.constraints_to_penalties
```

## The state vector and parameter object

Typically the unknowns of the system are present as a `Vector` of the appropriate length
in the numerical problem. The state vector can also be constructed manually without building
a problem.

The state vector `u` is the concatenation of the unknowns in the order of `unknowns(sys)`;
an array unknown `x` contributes `length(x)` consecutive entries, the generated functions
reconstruct it as a view into `u`, and symbolic indexing is unchanged (`prob[x]` returns
the whole array, `prob[x[i]]` one element).

```julia
@variables x(t)[1:3] y(t)
@parameters a[1:3]
eqs = [D(x[1]) ~ -a[1] * x[1], D(x[2]) ~ -a[2] * x[2], D(x[3]) ~ -a[3] * x[3] + y, D(y) ~ -y]
@named sys = System(eqs, t, [x, y], [a])
prob = ODEProblem(complete(sys), [x => ones(3), y => 1.0, a => [1.0, 2.0, 3.0]], (0.0, 1.0))
length(prob.u0) == 4
```

```@docs
ModelingToolkit.get_u0
ModelingToolkit.varmap_to_vars
```

The parameters of a split system are stored in a custom data structure called
`MTKParameters`. ModelingToolkit problem constructors use
[`SciMLBase.AutoDespecialize`](https://docs.sciml.ai/SciMLBase/stable/interfaces/Problems/)
by default. Solvers that support this policy wrap the parameters in
[`SciMLBase.DespecializedParameters`](https://docs.sciml.ai/SciMLBase/stable/interfaces/Problems/)
at solve time so compiled code
can be reused across parameter-buffer layouts. Explicit `AutoSpecialize` and
`FullSpecialize` problems retain their existing behavior. These objects should only be
interacted with through their defined public API.
SymbolicIndexingInterface.jl contains functionality useful for this purpose.

```@docs
MTKParameters
ModelingToolkit.get_p
```

The following functions are useful when working with `MTKParameters` objects, and especially
the `Tunables` portion. For more information about the "portions" of `MTKParameters`, refer
to the [`SciMLStructures.jl`](https://docs.sciml.ai/SciMLStructures/stable/) documentation.

```@docs
reorder_dimension_by_tunables!
reorder_dimension_by_tunables
```

## Initialization

```@docs
generate_initializesystem
InitializationProblem
ModelingToolkit.analyze_initialization_jacobian
ModelingToolkitBase.MissingGuessValue
```

## Linear analysis

```@docs
linearization_function
LinearizationProblem
ModelingToolkit.LinearizationOpPoint
linearize
ModelingToolkit.linearize_symbolic
```

There are also utilities for manipulating the results of these analyses in a symbolic context.

```@docs
ModelingToolkit.similarity_transform
ModelingToolkit.reorder_unknowns
```
