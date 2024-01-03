# Getting Started with ModelingToolkit.jl

This is an introductory tutorial for ModelingToolkit (MTK). We will demonstrate
the basics of the package by demonstrating how to define and simulate simple
Ordinary Differential Equation (ODE) systems.

## Installing ModelingToolkit

To install ModelingToolkit, use the Julia package manager. This can be done as follows:

```julia
using Pkg
Pkg.add("ModelingToolkit")
```

## Copy-Pastable Simplified Example

A much deeper tutorial with forcing functions and sparse Jacobians is below.
But if you want to just see some code and run, here's an example:

```@example first-mtkmodel
using ModelingToolkit

@variables t
D = Differential(t)

@mtkmodel FOL begin
    @parameters begin
        τ # parameters
    end
    @variables begin
        x(t) # dependent variables
    end
    @equations begin
        D(x) ~ (1 - x) / τ
    end
end

using DifferentialEquations: solve
@mtkbuild fol = FOL()
prob = ODEProblem(fol, [fol.x => 0.0], (0.0, 10.0), [fol.τ => 3.0])
sol = solve(prob)

using Plots
plot(sol)
```

Now let's start digging into MTK!

## Your very first ODE

Let us start with a minimal example. The system to be modelled is a
first-order lag element:

```math
\dot{x} = \frac{f(t) - x(t)}{\tau}
```

Here, ``t`` is the independent variable (time), ``x(t)`` is the (scalar) state
variable, ``f(t)`` is an external forcing function, and ``\tau`` is a
parameter.
In MTK, this system can be modelled as follows. For simplicity, we
first set the forcing function to a time-independent value ``1``. And the
independent variable ``t`` is automatically added by ``@mtkmodel``.

```@example ode2
using ModelingToolkit

@variables t
D = Differential(t)

@mtkmodel FOL begin
    @parameters begin
        τ # parameters
    end
    @variables begin
        x(t) # dependent variables
    end
    @equations begin
        D(x) ~ (1 - x) / τ
    end
end

@mtkbuild fol = FOL()
```

Note that equations in MTK use the tilde character (`~`) as equality sign.

`@mtkbuild` creates an instance of `FOL` named as `fol`.

After construction of the ODE, you can solve it using [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/):

```@example ode2
using DifferentialEquations
using Plots

prob = ODEProblem(fol, [fol.x => 0.0], (0.0, 10.0), [fol.τ => 3.0])
plot(solve(prob))
```

The initial state and the parameter values are specified using a mapping
from the actual symbolic elements to their values, represented as an array
of `Pair`s, which are constructed using the `=>` operator.

## Algebraic relations and structural simplification

You could separate the calculation of the right-hand side, by introducing an
intermediate variable `RHS`:

```@example ode2
using ModelingToolkit

@mtkmodel FOL begin
    @parameters begin
        τ # parameters
    end
    @variables begin
        x(t) # dependent variables
        RHS(t)
    end
    begin
        D = Differential(t)
    end
    @equations begin
        RHS ~ (1 - x) / τ
        D(x) ~ RHS
    end
end

@mtkbuild fol = FOL()
```

You can look at the equations by using the command `equations`:

```@example ode2
equations(fol)
```

Notice that there is only one equation in this system, `Differential(t)(x(t)) ~ RHS(t)`.
The other equation was removed from the system and was transformed into an `observed`
variable. Observed equations are variables which can be computed on-demand but are not
necessary for the solution of the system, and thus MTK tracks it separately. One can
check the observed equations via the `observed` function:

```@example ode2
observed(fol)
```

For more information on this process, see [Observables and Variable Elimination](@ref).

MTK still knows how to calculate them out of the information available
in a simulation result. The intermediate variable `RHS` therefore can be plotted
along with the state variable. Note that this has to be requested explicitly
like as follows:

```@example ode2
prob = ODEProblem(fol,
    [fol.x => 0.0],
    (0.0, 10.0),
    [fol.τ => 3.0])
sol = solve(prob)
plot(sol, vars = [fol.x, fol.RHS])
```

## Named Indexing of Solutions

Note that the indexing of the solution similarly works via the names, and so to get
the time series for `x`, one would do:

```@example ode2
sol[fol.x]
```

or to get the second value in the time series for `x`:

```@example ode2
sol[fol.x, 2]
```

Similarly, the time series for `RHS` can be retrieved using the same indexing:

```@example ode2
sol[fol.RHS]
```

## Specifying a time-variable forcing function

What if the forcing function (the “external input”) ``f(t)`` is not constant?
Obviously, one could use an explicit, symbolic function of time:

```@example ode2
@mtkmodel FOL begin
    @parameters begin
        τ # parameters
    end
    @variables begin
        x(t) # dependent variables
        f(t)
    end
    begin
        D = Differential(t)
    end
    @equations begin
        f ~ sin(t)
        D(x) ~ (f - x) / τ
    end
end

@named fol_variable_f = FOL()
```

But often this function might not be available in an explicit form.
Instead the function might be provided as time-series data.
MTK handles this situation by allowing us to “register” arbitrary Julia functions,
which are excluded from symbolic transformations, and thus used as-is.
So, you could, for example, interpolate given the time-series using
[DataInterpolations.jl](https://github.com/PumasAI/DataInterpolations.jl). Here,
we illustrate this option by a simple lookup ("zero-order hold") of a vector
of random values:

```@example ode2
value_vector = randn(10)
f_fun(t) = t >= 10 ? value_vector[end] : value_vector[Int(floor(t)) + 1]
@register_symbolic f_fun(t)

@mtkmodel FOLExternalFunction begin
    @parameters begin
        τ # parameters
    end
    @variables begin
        x(t) # dependent variables
        f(t)
    end
    @structural_parameters begin
        h = 1
    end
    begin
        D = Differential(t)
    end
    @equations begin
        f ~ f_fun(t)
        D(x) ~ (f - x) / τ
    end
end

@mtkbuild fol_external_f = FOLExternalFunction()
prob = ODEProblem(fol_external_f,
    [fol_external_f.x => 0.0],
    (0.0, 10.0),
    [fol_external_f.τ => 0.75])

sol = solve(prob)
plot(sol, vars = [fol_external_f.x, fol_external_f.f])
```

## Building component-based, hierarchical models

Working with simple one-equation systems is already fun, but composing more
complex systems from simple ones is even more fun. Best practice for such a
“modeling framework” could be to use factory functions for model components:

```@example ode2
function fol_factory(separate = false; name)
    @parameters τ
    @variables t x(t) f(t) RHS(t)

    eqs = separate ? [RHS ~ (f - x) / τ,
        D(x) ~ RHS] :
          D(x) ~ (f - x) / τ

    ODESystem(eqs; name)
end
```

Such a factory can then be used to instantiate the same component multiple times,
but allows for customization:

```@example ode2
@named fol_1 = fol_factory()
@named fol_2 = fol_factory(true) # has observable RHS
```

The `@named` macro rewrites `fol_2 = fol_factory(true)` into `fol_2 = fol_factory(true,:fol_2)`.
Now, these two components can be used as subsystems of a parent system, i.e.
one level higher in the model hierarchy. The connections between the components
again are just algebraic relations:

```@example ode2
connections = [fol_1.f ~ 1.5,
    fol_2.f ~ fol_1.x]

connected = compose(ODESystem(connections, name = :connected), fol_1, fol_2)
```

All equations, variables, and parameters are collected, but the structure of the
hierarchical model is still preserved. This means you can still get information about
`fol_1` by addressing it by `connected.fol_1`, or its parameter by
`connected.fol_1.τ`. Before simulation, we again eliminate the algebraic
variables and connection equations from the system using structural
simplification:

```@example ode2
connected_simp = structural_simplify(connected)
```

```@example ode2
full_equations(connected_simp)
```

As expected, only the two state-derivative equations remain,
as if you had manually eliminated as many variables as possible from the equations.
Some observed variables are not expanded unless `full_equations` is used.
As mentioned above, the hierarchical structure is preserved. So, the
initial state and the parameter values can be specified accordingly when
building the `ODEProblem`:

```@example ode2
u0 = [fol_1.x => -0.5,
    fol_2.x => 1.0]

p = [fol_1.τ => 2.0,
    fol_2.τ => 4.0]

prob = ODEProblem(connected_simp, u0, (0.0, 10.0), p)
plot(solve(prob))
```

More on this topic may be found in [Composing Models and Building Reusable Components](@ref acausal).

## Initial Guess

It is often a good idea to specify reasonable values for the initial state and the
parameters of a model component. Then, these do not have to be explicitly specified when constructing the `ODEProblem`.

```@example ode2
@mtkmodel UnitstepFOLFactory begin
    @parameters begin
        τ = 1.0
    end
    @variables begin
        x(t) = 0.0
    end
    @equations begin
        D(x) ~ (1 - x) / τ
    end
end
```

While defining the model `UnitstepFOLFactory`, an initial guess of 0.0 is assigned to `x(t)` and 1.0 to `τ`.
Additionally, these initial guesses can be modified while creating instances of `UnitstepFOLFactory` by passing arguments.

```@example ode2
@named fol = UnitstepFOLFactory(; x = 0.1)
sol = ODEProblem(fol, [], (0.0, 5.0), []) |> solve
```

In non-DSL definitions, one can pass `defaults` dictionary to set the initial guess of the symbolic variables.

```@example ode3
using ModelingToolkit

function UnitstepFOLFactory(; name)
    @parameters τ
    @variables t x(t)
    ODESystem(D(x) ~ (1 - x) / τ; name, defaults = Dict(x => 0.0, τ => 1.0))
end
```

Note that the defaults can be functions of the other variables, which is then
resolved at the time of the problem construction. Of course, the factory
function could accept additional arguments to optionally specify the initial
state or parameter values, etc.

## Symbolic and sparse derivatives

One advantage of a symbolic toolkit is that derivatives can be calculated
explicitly, and that the incidence matrix of partial derivatives (the
“sparsity pattern”) can also be explicitly derived. These two facts lead to a
substantial speedup of all model calculations, e.g. when simulating a model
over time using an ODE solver.

By default, analytical derivatives and sparse matrices, e.g. for the Jacobian, the
matrix of first partial derivatives, are not used. Let's benchmark this (`prob`
still is the problem using the `connected_simp` system above):

```@example ode2
using BenchmarkTools
@btime solve(prob, Rodas4());
nothing # hide
```

Now have MTK provide sparse, analytical derivatives to the solver. This has to
be specified during the construction of the `ODEProblem`:

```@example ode2
prob_an = ODEProblem(connected_simp, u0, (0.0, 10.0), p; jac = true)
@btime solve($prob_an, Rodas4());
nothing # hide
```

```@example ode2
prob_an = ODEProblem(connected_simp, u0, (0.0, 10.0), p; jac = true, sparse = true)
@btime solve($prob_an, Rodas4());
nothing # hide
```

The speedup is significant. For this small dense model (3 of 4 entries are
populated), using sparse matrices is counterproductive in terms of required
memory allocations. For large, hierarchically built models, which tend to be
sparse, speedup and the reduction of memory allocation can be expected to be
substantial. In addition, these problem builders allow for automatic parallelism
using the structural information. For more information, see the
[ODESystem](@ref ODESystem) page.

## Notes and pointers how to go on

Here are some notes that may be helpful during your initial steps with MTK:

  - The `@mtkmodel` macro is for high-level usage of MTK. However, in many cases you
    may need to programmatically generate `ODESystem`s. If that's the case, check out
    the [Programmatically Generating and Scripting ODESystems Tutorial](@ref programmatically).
  - Vector-valued parameters and variables are possible. A cleaner, more
    consistent treatment of these is still a work in progress, however. Once finished,
    this introductory tutorial will also cover this feature.

Where to go next?

  - Not sure how MTK relates to similar tools and packages? Read
    [Comparison of ModelingToolkit vs Equation-Based and Block Modeling Languages](@ref).
  - For a more detailed explanation of `@mtkmodel` checkout
    [Defining components with `@mtkmodel` and connectors with `@connectors`](@ref mtkmodel_connector)
  - Depending on what you want to do with MTK, have a look at some of the other
    **Symbolic Modeling Tutorials**.
  - If you want to automatically convert an existing function to a symbolic
    representation, you might go through the **ModelingToolkitize Tutorials**.
  - To learn more about the inner workings of MTK, consider the sections under
    **Basics** and **System Types**.
