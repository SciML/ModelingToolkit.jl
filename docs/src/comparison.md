# Comparison of Julia's ModelingToolkit vs SymPy for Symbolic Computation

ModelingToolkit.jl is a symbolic modeling language for Julia built in
Julia. Its goal is very different from Sympy: it was made to support
symbolic-numerics, the combination of symbolic computing with numerical
methods to allow for extreme performance computing that would not be
possible without modifying the model. Because of this, ModelingToolkit.jl
excels in many areas due to purposeful design decisions:

- Performance: ModelingToolkit.jl is built in Julia, whereas SymPy was
  built in Python. Thus the performance bar for ModelingToolkit.jl is
  much higher. ModelingToolkit.jl started because SymPy was far too
  slow and SymEngine was far too inflexible for the projects they were
  doing. Performance is key to ModelingToolkit.jl. If you find any
  performance issues, please file an issue.
- `build_function`: `lambdify` is "fine" for some people, but if you're building
  a super fast MPI-enabled Julia/C/Fortran simulation code, having a
  function that hits the Python interpreter is less than optimal. By
  default, `build_function` builds fast JIT-compiled functions due
  to being in Julia. However, it has support for things like static
  arrays, non-allocating functions via mutation, fast functions on
  sparse matrices and arrays of arrays, etc.: all core details of
  doing high performance computing.
- Parallelism: ModelingToolkit.jl has pervasive parallelism. The
  symbolic simplification via [SymbolicUtils.jl](https://github.com/JuliaSymbolics/SymbolicUtils.jl)
  has built-in parallelism, ModelingToolkit.jl builds functions that
  parallelizes across threads. ModelingToolkit.jl is compatible with GPU libraries like CUDA.jl.
- Scientific Machine Learning (SciML): ModelingToolkit.jl is made to synergize
  with the high performance Julia SciML ecosystem in many ways. At a
  base level, all expressions and built functions are compatible with
  automatic differentiation like ForwardDiff.jl and Zygote.jl, meaning
  that it can be used in and with neural networks. Tools like
  [DataDrivenDiffEq.jl](https://datadriven.sciml.ai/dev/) can reconstruct
  symbolic expressions from neural networks and data while
  [NeuralNetDiffEq.jl](https://github.com/SciML/NeuralNetDiffEq.jl)
  can automatically solve partial differential equations from symbolic
  descriptions using physics-informed neural networks.
- Primitives for high-performance numerics. Features like `ODESystem`
  can be used to easily generate automatically parallelized ODE solver
  code with sparse Jacobians and all of the pieces required to get
  the most optimal solves. Support for differential-algebraic equations,
  chemical reaction networks, and generation of code for nonlinear
  optimization tools makes ModelingToolkit.jl a tool for, well,
  building, generating, and analyzing models.
- Deep integration with the Julia ecosystem: ModelingToolkit.jl's integration
  with neural networks is not the only thing that's deep. ModelingToolkit.jl
  is built with the same philosophy as other SciML packages, eschewing
  "monorepos" for a distributed development approach that ties together
  the work of many developers. The differentiation parts utilize tools
  from automatic differentiation libraries, all linear algebra functionality
  comes from tracing Julia Base itself, symbolic rewriting (simplification
  and substitution) comes from
  [SymbolicUtils.jl](https://github.com/JuliaSymbolics/SymbolicUtils.jl),
  parallelism comes from Julia Base libraries and Dagger.jl, and etc.
  The list keeps going. All told, by design ModelingToolkit.jl's development
  moves fast because it's effectively using the work of hundreds of
  Julia developers, allowing it to grow fast.
