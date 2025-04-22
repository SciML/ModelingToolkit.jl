# [Programmatically Generating and Scripting ODESystems](@id programmatically)

In the following tutorial, we will discuss how to programmatically generate `ODESystem`s.
This is useful for functions that generate `ODESystem`s, for example
when you implement a reader that parses some file format, such as SBML, to generate an `ODESystem`.
It is also useful for functions that transform an `ODESystem`, for example
when you write a function that log-transforms a variable in an `ODESystem`.

## The Representation of a ModelingToolkit System

ModelingToolkit is built on [Symbolics.jl](https://symbolics.juliasymbolics.org/dev/),
a symbolic Computer Algebra System (CAS) developed in Julia. As such, all CAS functionality
is also available to be used on ModelingToolkit systems, such as symbolic differentiation, Groebner basis
calculations, and whatever else you can think of. Under the hood, all ModelingToolkit
variables and expressions are Symbolics.jl variables and expressions. Thus when scripting
a ModelingToolkit system, one simply needs to generate Symbolics.jl variables and equations
as demonstrated in the Symbolics.jl documentation. This looks like:

```@example scripting
using ModelingToolkit # reexports Symbolics
@variables t x(t) y(t) # Define variables
D = Differential(t)
eqs = [D(x) ~ y
       D(y) ~ x] # Define an array of equations
```

However, ModelingToolkit has many higher-level features which will make scripting ModelingToolkit systems more convenient.
For example, as shown in the next section, defining your own independent variables and differentials is rarely needed.

## The Non-DSL (non-`@mtkmodel`) Way of Defining an ODESystem

Using `@mtkmodel`, like in the [getting started tutorial](@ref getting_started),
is the preferred way of defining ODEs with MTK.
However generating the contents of a `@mtkmodel` programmatically can be tedious.
Let us look at how we can define the same system without `@mtkmodel`.

```@example scripting
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
@variables x(t) = 0.0  # independent and dependent variables
@parameters τ = 3.0       # parameters
@constants h = 1    # constants
eqs = [D(x) ~ (h - x) / τ] # create an array of equations

# your first ODE, consisting of a single equation, indicated by ~
@named model = ODESystem(eqs, t)

# Perform the standard transformations and mark the model complete
# Note: Complete models cannot be subsystems of other models!
fol = structural_simplify(model)
prob = ODEProblem(fol, [], (0.0, 10.0), [])
using OrdinaryDiffEq
sol = solve(prob)

using Plots
plot(sol)
```

As you can see, generating an ODESystem is as simple as creating an array of equations
and passing it to the `ODESystem` constructor.

`@named` automatically gives a name to the `ODESystem`, and is shorthand for

```@example scripting
fol_model = ODESystem(eqs, t; name = :fol_model) # @named fol_model = ODESystem(eqs, t)
```

Thus, if we had read a name from a file and wish to populate an `ODESystem` with said name, we could do:

```@example scripting
namesym = :name_from_file
fol_model = ODESystem(eqs, t; name = namesym)
```

## Warning About Mutation

Be advsied that it's never a good idea to mutate an `ODESystem`, or any `AbstractSystem`.
