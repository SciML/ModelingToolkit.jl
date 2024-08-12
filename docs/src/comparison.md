# Comparison of ModelingToolkit vs Equation-Based and Block Modeling Languages

## Comparison Against Modelica

  - Both Modelica and ModelingToolkit.jl are acausal modeling languages.
  - Modelica is a language with many different implementations, such as
    [Dymola](https://www.3ds.com/products/catia/dymola/) and
    [OpenModelica](https://openmodelica.org/), which have differing levels of
    performance and can give different results on the same model. Many of the
    commonly used Modelica compilers are not open-source. ModelingToolkit.jl
    is a language with a single canonical open-source implementation.
  - All current Modelica compiler implementations are fixed and not extendable
    by the users from the Modelica language itself. For example, the Dymola
    compiler [shares its symbolic processing pipeline](https://www.claytex.com/tech-blog/model-translation-and-symbolic-manipulation/),
    which is roughly equivalent to the `dae_index_lowering` and `structural_simplify`
    of ModelingToolkit.jl. ModelingToolkit.jl is an open and hackable transformation
    system which allows users to add new non-standard transformations and control
    the order of application.
  - Modelica is a declarative programming language. ModelingToolkit.jl is a
    declarative symbolic modeling language used from within the Julia programming
    language. Its programming language semantics, such as loop constructs and
    conditionals, can be used to more easily generate models.
  - Modelica is an object-oriented single dispatch language. ModelingToolkit.jl,
    built on Julia, uses multiple dispatch extensively to simplify code.
  - Many Modelica compilers supply a GUI. ModelingToolkit.jl does not.
  - Modelica can be used to simulate ODE and DAE systems. ModelingToolkit.jl
    has a much more expansive set of system types, including nonlinear systems,
    SDEs, PDEs, and more.

## Comparison Against Simulink

  - Simulink is a causal modeling environment, whereas ModelingToolkit.jl is an
    acausal modeling environment. For an overview of the differences, consult
    academic reviews such as [this one](https://arxiv.org/abs/1909.00484). In this
    sense, ModelingToolkit.jl is more similar to the Simscape sub-environment.
  - Simulink is used from MATLAB while ModelingToolkit.jl is used from Julia.
    Thus any user-defined functions have the performance of their host language.
    For information on the performance differences between Julia and MATLAB,
    consult [open benchmarks](https://julialang.org/benchmarks/), which demonstrate
    Julia as an order of magnitude or more faster in many cases due to its JIT
    compilation.
  - Simulink uses the MATLAB differential equation solvers, while ModelingToolkit.jl
    uses [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/). For a systematic
    comparison between the solvers, consult
    [open benchmarks](https://docs.sciml.ai/SciMLBenchmarksOutput/stable/),
    which demonstrate two orders of magnitude performance advantage for the native
    Julia solvers across many benchmark problems.
  - Simulink comes with a Graphical User Interface (GUI), ModelingToolkit.jl
    does not.
  - Simulink is a proprietary software, meaning users cannot actively modify or
    extend the software. ModelingToolkit.jl is built in Julia and used in Julia,
    where users can actively extend and modify the software interactively in the
    REPL and contribute to its open-source repositories.
  - Simulink covers ODE and DAE systems. ModelingToolkit.jl has a much more
    expansive set of system types, including SDEs, PDEs, optimization problems,
    and more.

## Comparison Against CASADI

  - CASADI is written in C++ but used from Python/MATLAB, meaning that it cannot be
    directly extended by users unless they are using the C++ interface and run a
    local build of CASADI. ModelingToolkit.jl is both written and used from
    Julia, meaning that users can easily extend the library on the fly, even
    interactively in the REPL.
  - CASADI includes limited support for Computer Algebra System (CAS) functionality,
    while ModelingToolkit.jl is built on the full
    [Symbolics.jl](https://docs.sciml.ai/Symbolics/stable/) CAS.
  - CASADI supports DAE and ODE problems via SUNDIALS IDAS and CVODES. ModelingToolkit.jl
    supports DAE and ODE problems via [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/),
    of which Sundials.jl is <1% of the total available solvers and is outperformed
    by the native Julia solvers on the vast majority of the benchmark equations.
    In addition, the DifferentialEquations.jl interface is confederated, meaning
    that any user can dynamically extend the system to add new solvers to the
    interface by defining new dispatches of solve.
  - CASADI's DAEBuilder does not implement efficiency transformations like tearing,
    which are standard in the ModelingToolkit.jl transformation pipeline.
  - CASADI supports special functionality for quadratic programming problems, while
    ModelingToolkit only provides nonlinear programming via `OptimizationSystem`.
  - ModelingToolkit.jl integrates with its host language Julia, so Julia code
    can be automatically converted into ModelingToolkit expressions. Users of
    CASADI must explicitly create CASADI expressions.

## Comparison Against Modia.jl

  - Modia.jl uses Julia's expression objects for representing its equations.
    ModelingToolkit.jl uses [Symbolics.jl](https://docs.sciml.ai/Symbolics/stable/),
    and thus the Julia expressions follow Julia semantics and can be manipulated
    using a computer algebra system (CAS).
  - Modia's compilation pipeline is similar to the
    [Dymola symbolic processing pipeline](https://www.claytex.com/tech-blog/model-translation-and-symbolic-manipulation/)
    with some improvements. ModelingToolkit.jl has an open transformation pipeline
    that allows for users to extend and reorder transformation passes, where
    `structural_simplify` is an adaptation of the Modia.jl-improved alias elimination
    and tearing algorithms.
  - Both Modia and ModelingToolkit generate `DAEProblem` and `ODEProblem` forms for
    solving with [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/).
  - ModelingToolkit.jl integrates with its host language Julia, so Julia code
    can be automatically converted into ModelingToolkit expressions. Users of
    Modia must explicitly create Modia expressions.
  - Modia covers DAE systems. ModelingToolkit.jl has a much more
    expansive set of system types, including SDEs, PDEs, optimization problems,
    and more.

## Comparison Against Causal.jl

  - Causal.jl is a causal modeling environment, whereas ModelingToolkit.jl is an
    acausal modeling environment. For an overview of the differences, consult
    academic reviews such as [this one](https://arxiv.org/abs/1909.00484).
  - Both ModelingToolkit.jl and Causal.jl use [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/)
    as the backend solver library.
  - Causal.jl lets one add arbitrary equation systems to a given node, and allow
    the output to effect the next node. This means an SDE may drive an ODE. These
    two portions are solved with different solver methods in tandem. In
    ModelingToolkit.jl, such connections promote the whole system to an SDE. This
    results in better accuracy and stability, though in some cases it can be
    less performant.
  - Causal.jl, similar to Simulink, breaks algebraic loops via inexact heuristics.
    ModelingToolkit.jl treats algebraic loops exactly through algebraic equations
    in the generated model.
