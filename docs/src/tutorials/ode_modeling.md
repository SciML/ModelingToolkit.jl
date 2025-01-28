# [Getting Started with ModelingToolkit.jl](@id getting_started)

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
But if you want to just see some code and run it, here's an example:

```@example first-mtkmodel
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@mtkmodel FOL begin
    @parameters begin
        τ = 3.0 # parameters
    end
    @variables begin
        x(t) = 0.0 # dependent variables
    end
    @equations begin
        D(x) ~ (1 - x) / τ
    end
end

using OrdinaryDiffEq
@mtkbuild fol = FOL()
prob = ODEProblem(fol, [], (0.0, 10.0), [])
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

Here, ``t`` is the independent variable (time), ``x(t)`` is the (scalar) unknown
variable, ``f(t)`` is an external forcing function, and ``\tau`` is a
parameter.
In MTK, this system can be modelled as follows. For simplicity, we
first set the forcing function to a time-independent value ``1``. And the
independent variable ``t`` is automatically added by `@mtkmodel`.

```@example ode2
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@mtkmodel FOL begin
    @parameters begin
        τ = 3.0 # parameters and their values
    end
    @variables begin
        x(t) = 0.0 # dependent variables and their initial conditions
    end
    @equations begin
        D(x) ~ (1 - x) / τ
    end
end

@mtkbuild fol = FOL()
```

Note that equations in MTK use the tilde character (`~`) as equality sign.

`@mtkbuild` creates an instance of `FOL` named as `fol`.

After construction of the ODE, you can solve it using [OrdinaryDiffEq.jl](https://docs.sciml.ai/DiffEqDocs/stable/):

```@example ode2
using OrdinaryDiffEq
using Plots

prob = ODEProblem(fol, [], (0.0, 10.0), [])
plot(solve(prob))
```

The parameter values are determined using the right hand side of the expressions in the `@parameters` block,
and similarly initial conditions are determined using the right hand side of the expressions in the `@variables` block.

## Using different values for parameters and initial conditions

If you want to simulate the same model,
but with different values for the parameters and initial conditions than the default values,
you likely do not want to write an entirely new `@mtkmodel`.
ModelingToolkit supports overwriting the default values:

```@example ode2
@mtkbuild fol_different_values = FOL(; τ = 1 / 3, x = 0.5)
prob = ODEProblem(fol_different_values, [], (0.0, 10.0), [])
plot(solve(prob))
```

Alternatively, this overwriting could also have occurred at the `ODEProblem` level.

```@example ode2
prob = ODEProblem(fol, [fol.τ => 1 / 3], (0.0, 10.0), [fol.x => 0.5])
plot(solve(prob))
```

Here, the second argument of `ODEProblem` is an array of `Pairs`.
The left hand side of each Pair is the parameter you want to overwrite,
and the right hand side is the value to overwrite it with.
Similarly, the initial conditions are overwritten in the fourth argument.
One important difference with the previous method is
that the parameter has to be referred to as `fol.τ` instead of just `τ`.

## Algebraic relations and structural simplification

You could separate the calculation of the right-hand side, by introducing an
intermediate variable `RHS`:

```@example ode2
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@mtkmodel FOL begin
    @parameters begin
        τ = 3.0 # parameters and their values
    end
    @variables begin
        x(t) = 0.0 # dependent variables and their initial conditions
        RHS(t)
    end
    @equations begin
        RHS ~ (1 - x) / τ
        D(x) ~ RHS
    end
end

@mtkbuild fol = FOL()
```

If you copy this block of code to your REPL, you will not see the above LaTeX equations.
Instead, you can look at the equations by using the `equations` function:

```@example ode2
equations(fol)
```

Notice that there is only one equation in this system, `Differential(t)(x(t)) ~ RHS(t)`.
The other equation was removed from the system and was transformed into an `observed`
variable. Observed equations are variables that can be computed on-demand but are not
necessary for the solution of the system, and thus MTK tracks them separately.
For this reason, we also did not need to specify an initial condition for `RHS`.
You can check the observed equations via the `observed` function:

```@example ode2
observed(fol)
```

For more information on this process, see [Observables and Variable Elimination](@ref).

MTK still knows how to calculate them out of the information available
in a simulation result. The intermediate variable `RHS` therefore can be plotted
along with the unknown variable. Note that this has to be requested explicitly:

```@example ode2
prob = ODEProblem(fol, [], (0.0, 10.0), [])
sol = solve(prob)
plot(sol, idxs = [fol.x, fol.RHS])
```

## Named Indexing of Solutions

Note that the indexing of the solution also works via the symbol, and so to get
the time series for `x`, you would do:

```@example ode2
sol[fol.x]
```

or to get the second value in the time series for `x`:

```@example ode2
sol[fol.x, 2]
```

Similarly, the time series for `RHS` can be retrieved using the same symbolic indexing:

```@example ode2
sol[fol.RHS]
```

## Specifying a time-variable forcing function

What if the forcing function (the “external input”) ``f(t)`` is not constant?
Obviously, one could use an explicit, symbolic function of time:

```@example ode2
@mtkmodel FOL begin
    @parameters begin
        τ = 3.0 # parameters and their values
    end
    @variables begin
        x(t) = 0.0 # dependent variables and their initial conditions
        f(t)
    end
    @equations begin
        f ~ sin(t)
        D(x) ~ (f - x) / τ
    end
end

@mtkbuild fol_variable_f = FOL()
```

However, this function might not be available in an explicit form.
Instead, the function might be provided as time-series data.
MTK handles this situation by allowing us to “register” arbitrary Julia functions,
which are excluded from symbolic transformations and thus used as-is.
For example, you could interpolate given the time-series using
[DataInterpolations.jl](https://github.com/SciML/DataInterpolations.jl). Here,
we illustrate this option with a simple lookup ("zero-order hold") of a vector
of random values:

```@example ode2
value_vector = randn(10)
f_fun(t) = t >= 10 ? value_vector[end] : value_vector[Int(floor(t)) + 1]
@register_symbolic f_fun(t)

@mtkmodel FOLExternalFunction begin
    @parameters begin
        τ = 0.75 # parameters and their values
    end
    @variables begin
        x(t) = 0.0 # dependent variables and their initial conditions
        f(t)
    end
    @equations begin
        f ~ f_fun(t)
        D(x) ~ (f - x) / τ
    end
end

@mtkbuild fol_external_f = FOLExternalFunction()
```

```@example ode2
prob = ODEProblem(fol_external_f, [], (0.0, 10.0), [])
sol = solve(prob)
plot(sol, idxs = [fol_external_f.x, fol_external_f.f])
```

## Building component-based, hierarchical models

Working with simple one-equation systems is already fun, but composing more
complex systems from simple ones is even more fun. The best practice for such a
“modeling framework” is to use the `@components` block in the `@mtkmodel` macro:

```@example ode2
@mtkmodel FOLUnconnectedFunction begin
    @parameters begin
        τ # parameters
    end
    @variables begin
        x(t) # dependent variables
        f(t)
        RHS(t)
    end
    @equations begin
        RHS ~ f
        D(x) ~ (RHS - x) / τ
    end
end
@mtkmodel FOLConnected begin
    @components begin
        fol_1 = FOLUnconnectedFunction(; τ = 2.0, x = -0.5)
        fol_2 = FOLUnconnectedFunction(; τ = 4.0, x = 1.0)
    end
    @equations begin
        fol_1.f ~ 1.5
        fol_2.f ~ fol_1.x
    end
end
@mtkbuild connected = FOLConnected()
```

Here the total model consists of two of the same submodels (components),
but with a different input function, parameter values and initial conditions.
The first model has a constant input, and the second model uses the state `x` of the first system as an input.
To avoid having to type the same differential equation multiple times,
we define the submodel in a separate `@mtkmodel`.
We then reuse this submodel twice in the total model `@components` block.
The inputs of two submodels then still have to be specified in the `@equations` block.

All equations, variables, and parameters are collected, but the structure of the
hierarchical model is still preserved. This means you can still get information about
`fol_1` by addressing it by `connected.fol_1`, or its parameter by
`connected.fol_1.τ`.

As expected, only the two equations with the derivatives of unknowns remain,
as if you had manually eliminated as many variables as possible from the equations.
Some observed variables are not expanded unless `full_equations` is used.
As mentioned above, the hierarchical structure is preserved. So, the
initial unknown and the parameter values can be specified accordingly when
building the `ODEProblem`:

```@example ode2
prob = ODEProblem(connected, [], (0.0, 10.0), [])
plot(solve(prob))
```

More on this topic may be found in [Composing Models and Building Reusable Components](@ref acausal).

## Symbolic and sparse derivatives

One advantage of a symbolic toolkit is that derivatives can be calculated
explicitly, and that the incidence matrix of partial derivatives (the
“sparsity pattern”) can also be explicitly derived. These two facts lead to a
substantial speedup of all model calculations, e.g. when simulating a model
over time using an ODE solver.

By default, analytical derivatives and sparse matrices, e.g. for the Jacobian, the
matrix of first partial derivatives, are not used. Let's benchmark this (`prob`
still is the problem using the `connected` system above):

```@example ode2
using BenchmarkTools
@btime solve(prob, Rodas4());
nothing # hide
```

Now have MTK provide sparse, analytical derivatives to the solver. This has to
be specified during the construction of the `ODEProblem`:

```@example ode2
prob_an = ODEProblem(connected, [], (0.0, 10.0), []; jac = true)
@btime solve(prob_an, Rodas4());
nothing # hide
```

```@example ode2
prob_sparse = ODEProblem(connected, [], (0.0, 10.0), []; jac = true, sparse = true)
@btime solve(prob_sparse, Rodas4());
nothing # hide
```

The speedup using the analytical Jacobian is significant.
For this small dense model (3 of 4 entries populated),
using sparse matrices is counterproductive in terms of required
memory allocations. For large, hierarchically built models, which tend to be
sparse, speedup and the reduction of memory allocation can also be expected to be
substantial. In addition, these problem builders allow for automatic parallelism by
exploiting the structural information. For more information, see the
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
    [Defining components with `@mtkmodel` and connectors with `@connectors`](@ref mtk_language)
  - Depending on what you want to do with MTK, have a look at some of the other
    **Symbolic Modeling Tutorials**.
  - If you want to automatically convert an existing function to a symbolic
    representation, you might go through the **ModelingToolkitize Tutorials**.
  - To learn more about the inner workings of MTK, consider the sections under
    **Basics** and **System Types**.
