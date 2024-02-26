# [Programmatically Generating and Scripting ODESystems](@id programmatically)

In the following tutorial we will discuss how to programmatically generate `ODESystem`s.
This is for cases where one is writing functions that generating `ODESystem`s, for example
if implementing a reader which parses some file format to generate an `ODESystem` (for example,
SBML), or for writing functions that transform an `ODESystem` (for example, if you write a
function that log-transforms a variable in an `ODESystem`).

## The Representation of a ModelingToolkit System

ModelingToolkit is built on [Symbolics.jl](https://symbolics.juliasymbolics.org/dev/),
a symbolic Computer Algebra System (CAS) developed in Julia. As such, all CAS functionality
is available on ModelingToolkit systems, such as symbolic differentiation, Groebner basis
calculations, and whatever else you can think of. Under the hood, all ModelingToolkit
variables and expressions are Symbolics.jl variables and expressions. Thus when scripting
a ModelingToolkit system, one simply needs to generate Symbolics.jl variables and equations
as demonstrated in the Symbolics.jl documentation. This looks like:

```@example scripting
using Symbolics
using ModelingToolkit: t_nounits as t, D_nounits as D

@variables x(t) y(t) # Define variables
eqs = [D(x) ~ y
       D(y) ~ x] # Define an array of equations
```

## The Non-DSL (non-`@mtkmodel`) Way of Defining an ODESystem

Using `@mtkmodel` is the preferred way of defining ODEs with MTK. However, let us
look at how we can define the same system without `@mtkmodel`. This is useful for
defining PDESystem etc.

```@example scripting
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@variables x(t)   # independent and dependent variables
@parameters τ       # parameters
@constants h = 1    # constants
eqs = [D(x) ~ (h - x) / τ] # create an array of equations

# your first ODE, consisting of a single equation, indicated by ~
@named model = ODESystem(eqs, t)

# Perform the standard transformations and mark the model complete
# Note: Complete models cannot be subsystems of other models!
fol_model = structural_simplify(model)
```

As you can see, generating an ODESystem is as simple as creating the array of equations
and passing it to the `ODESystem` constructor.

## Understanding the Difference Between the Julia Variable and the Symbolic Variable

In the most basic usage of ModelingToolkit and Symbolics, the name of the Julia variable
and the symbolic variable are the same. For example, when we do:

```@example scripting
@variables a
```

the name of the symbolic variable is `a` and same with the Julia variable. However, we can
de-couple these by setting `a` to a new symbolic variable, for example:

```@example scripting
b = only(@variables(a))
```

Now the Julia variable `b` refers to the variable named `a`. However, the downside of this current
approach is that it requires that the user writing the script knows the name `a` that they want to
place to the variable. But what if for example we needed to get the variable's name from a file?

To do this, one can interpolate a symbol into the `@variables` macro using `$`. For example:

```@example scripting
a = :c
b = only(@variables($a))
```

In this example, `@variables($a)` created a variable named `c`, and set this variable to `b`.

Variables are not the only thing with names. For example, when you build a system, it knows its name
that name is used in the namespacing. In the standard usage, again the Julia variable and the
symbolic name are made the same via:

```@example scripting
@named fol_model = ODESystem(eqs, t)
```

However, one can decouple these two properties by noting that `@named` is simply shorthand for the
following:

```@example scripting
fol_model = ODESystem(eqs, t; name = :fol_model)
```

Thus if we had read a name from a file and wish to populate an `ODESystem` with said name, we could do:

```@example scripting
namesym = :name_from_file
fol_model = ODESystem(eqs, t; name = namesym)
```

## Warning About Mutation

Be advsied that it's never a good idea to mutate an `ODESystem`, or any `AbstractSystem`.
