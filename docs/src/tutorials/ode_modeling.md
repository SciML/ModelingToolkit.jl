# Getting Started with ModelingToolkit.jl

This is an introductory tutorial for ModelingToolkit (MTK).
Some examples of Ordinary Differential Equations (ODE) are used to
illustrate the basic user-facing functionality.

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
    @structural_parameters begin
        h = 1
    end
    @equations begin
        D(x) ~ (h - x) / τ
    end
end

using DifferentialEquations: solve

@named fol = FOL()
fol = complete(fol)

prob = ODEProblem(fol, [fol.x => 0.0], (0.0, 10.0), [fol.τ => 3.0])
# parameter `τ` can be assigned a value, but structural parameter `h` cannot'.
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
first set the forcing function to a time-independent value ``h``. And the
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
    @structural_parameters begin
        h = 1
    end
    @equations begin
        D(x) ~ (h - x) / τ
    end
end

@named fol_incomplete = FOL()
fol = complete(fol_incomplete)
```

Note that equations in MTK use the tilde character (`~`) as equality sign.

`@named` creates an instance of `FOL` named as `fol`. Before creating an
ODEProblem with `fol` run `complete`. Once the system is complete, it will no
longer namespace its subsystems or variables. This is necessary to correctly pass
the intial values of states and parameters to the ODEProblem.

```julia
julia> fol_incomplete.x
fol_incomplete₊x(t)

julia> fol.x
x(t)
```

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

## Non-DSL way of defining an ODESystem

Using `@mtkmodel` is the preferred way of defining ODEs with MTK. However, let us
look at how we can define the same system without `@mtkmodel`. This is useful for
defining PDESystem etc.

```@example first-mtkmodel
@variables t x(t)   # independent and dependent variables
@parameters τ       # parameters
@constants h = 1    # constants
D = Differential(t) # define an operator for the differentiation w.r.t. time

# your first ODE, consisting of a single equation, indicated by ~
@named fol_model = ODESystem(D(x) ~ (h - x) / τ)
```

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
    @structural_parameters begin
        h = 1
    end
    begin
        D = Differential(t)
    end
    @equations begin
        RHS ~ (h - x) / τ
        D(x) ~ RHS
    end
end

@named fol_separate = FOL()
```

To directly solve this system, you would have to create a Differential-Algebraic
Equation (DAE) problem, since besides the differential equation, there is an
additional algebraic equation now. However, this DAE system can obviously be
transformed into the single ODE we used in the first example above. MTK achieves
this by structural simplification:

```@example ode2
fol_simplified = structural_simplify(complete(fol_separate))
equations(fol_simplified)
```

```@example ode2
equations(fol_simplified) == equations(fol)
```

You can extract the equations from a system using `equations` (and, in the same
way, `states` and `parameters`). The simplified equation is exactly the same
as the original one, so the simulation performance will also be the same.
However, there is one difference. MTK does keep track of the eliminated
algebraic variables as "observables" (see
[Observables and Variable Elimination](@ref)).
That means, MTK still knows how to calculate them out of the information available
in a simulation result. The intermediate variable `RHS` therefore can be plotted
along with the state variable. Note that this has to be requested explicitly,
through:

```@example ode2
prob = ODEProblem(fol_simplified,
    [fol_simplified.x => 0.0],
    (0.0, 10.0),
    [fol_simplified.τ => 3.0])
sol = solve(prob)
plot(sol, vars = [fol_simplified.x, fol_simplified.RHS])
```

By default, `structural_simplify` also replaces symbolic `constants` with
their default values. This allows additional simplifications not possible
when using `parameters` (e.g., solution of linear equations by dividing out
the constant's value, which cannot be done for parameters, since they may
be zero).

Note that the indexing of the solution similarly works via the names, and so
`sol[x]` gives the time-series for `x`, `sol[x,2:10]` gives the 2nd through 10th
values of `x` matching `sol.t`, etc. Note that this works even for variables
which have been eliminated, and thus `sol[RHS]` retrieves the values of `RHS`.

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
    @structural_parameters begin
        h = 1
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

@named fol_external_f = FOLExternalFunction()
fol_external_f = complete(fol_external_f)
prob = ODEProblem(structural_simplify(fol_external_f),
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

## Inital Guess

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
Additionaly, these initial guesses can be modified while creating instances of `UnitstepFOLFactory` by passing arguements.

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

  - Sometimes, the symbolic engine within MTK cannot correctly identify the
    independent variable (e.g. time) out of all variables. In such a case, you
    usually get an error that some variable(s) is "missing from variable map". In
    most cases, it is then sufficient to specify the independent variable as second
    argument to `ODESystem`, e.g. `ODESystem(eqs, t)`.
  - A completely macro-free usage of MTK is possible and is discussed in a
    separate tutorial. This is for package developers, since the macros are only
    essential for automatic symbolic naming for modelers.
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
