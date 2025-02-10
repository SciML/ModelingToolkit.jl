# [Symbolic Metadata](@id symbolic_metadata)

It is possible to add metadata to symbolic variables, the metadata will be displayed when calling help on a variable.

The following information can be added (note, it's possible to extend this to user-defined metadata as well)

## Variable descriptions

Descriptive strings can be attached to variables using the `[description = "descriptive string"]` syntax:

```@example metadata
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
@variables u [description = "This is my input"]
getdescription(u)
```

When variables with descriptions are present in systems, they will be printed when the system is shown in the terminal:

```@example metadata
@variables u(t) [description = "A short description of u"]
@parameters p [description = "A description of p"]
@named sys = ODESystem([u ~ p], t)
show(stdout, "text/plain", sys) # hide
```

Calling help on the variable `u` displays the description, alongside other metadata:

```
help?> u

  A variable of type Symbolics.Num (Num wraps anything in a type that is a subtype of Real)

  Metadata
  ≡≡≡≡≡≡≡≡≡≡

  ModelingToolkit.VariableDescription: This is my input

  Symbolics.VariableSource: (:variables, :u)
```

## Connect

Variables in connectors can have `connect` metadata which describes the type of connections.

`Flow` is used for variables that represent physical quantities that "flow" ex:
current in a resistor. These variables sum up to zero in connections.

`Stream` can be specified for variables that flow bi-directionally.

```@example connect
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@variables i(t) [connect = Flow]
@variables k(t) [connect = Stream]
hasconnect(i)
```

```@example connect
getconnect(k)
```

## Input or output

Designate a variable as either an input or an output using the following

```@example metadata
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@variables u [input = true]
isinput(u)
```

```@example metadata
@variables y [output = true]
isoutput(y)
```

## Bounds

Bounds are useful when parameters are to be optimized, or to express intervals of uncertainty.

```@example metadata
@variables u [bounds = (-1, 1)]
hasbounds(u)
```

```@example metadata
getbounds(u)
```

Bounds can also be specified for array variables. A scalar array bound is applied to each
element of the array. A bound may also be specified as an array, in which case the size of
the array must match the size of the symbolic variable.

```@example metadata
@variables x[1:2, 1:2] [bounds = (-1, 1)]
hasbounds(x)
```

```@example metadata
getbounds(x)
```

```@example metadata
getbounds(x[1, 1])
```

```@example metadata
getbounds(x[1:2, 1])
```

```@example metadata
@variables x[1:2] [bounds = (-Inf, [1.0, Inf])]
hasbounds(x)
```

```@example metadata
getbounds(x)
```

```@example metadata
getbounds(x[2])
```

```@example metadata
hasbounds(x[2])
```

## Guess

Specify an initial guess for custom initial conditions of an `ODESystem`.

```@example metadata
@variables u [guess = 1]
hasguess(u)
```

```@example metadata
getguess(u)
```

## Mark input as a disturbance

Indicate that an input is not available for control, i.e., it's a disturbance input.

```@example metadata
@variables u [input = true, disturbance = true]
isdisturbance(u)
```

## Mark parameter as tunable

Indicate that a parameter can be automatically tuned by parameter optimization or automatic control tuning apps.

```@example metadata
@parameters Kp [tunable = true]
istunable(Kp)
```

## Probability distributions

A probability distribution may be associated with a parameter to indicate either
uncertainty about its value, or as a prior distribution for Bayesian optimization.

```julia
using Distributions
d = Normal(10, 1)
@parameters m [dist = d]
hasdist(m)
```

```julia
getdist(m)
```

## Irreducible

A variable can be marked `irreducible` to prevent it from being moved to an
`observed` state. This forces the variable to be computed during solving so that
it can be accessed in [callbacks](@ref events)

```@example metadata
@variables important_value [irreducible = true]
isirreducible(important_value)
```

## State Priority

When a model is structurally simplified, the algorithm will try to ensure that the variables with higher state priority become states of the system. A variable's state priority is a number set using the `state_priority` metadata.

```@example metadata
@variables important_dof [state_priority = 10] unimportant_dof [state_priority = -2]
state_priority(important_dof)
```

## Units

Units for variables can be designated using symbolic metadata. For more information, please see the [model validation and units](@ref units) section of the docs. Note that `getunit` is not equivalent to `get_unit` - the former is a metadata getter for individual variables (and is provided so the same interface function for `unit` exists like other metadata), while the latter is used to handle more general symbolic expressions.

```@example metadata
@variables speed [unit = u"m/s"]
hasunit(speed)
```

```@example metadata
getunit(speed)
```

## Miscellaneous metadata

User-defined metadata can be added using the `misc` metadata. This can be queried
using the `hasmisc` and `getmisc` functions.

```@example metadata
@variables u [misc = :conserved_parameter] y [misc = [2, 4, 6]]
hasmisc(u)
```

```@example metadata
getmisc(y)
```

## Additional functions

For systems that contain parameters with metadata like described above, have some additional functions defined for convenience.
In the example below, we define a system with tunable parameters and extract bounds vectors

```@example metadata
@variables x(t)=0 u(t)=0 [input = true] y(t)=0 [output = true]
@parameters T [tunable = true, bounds = (0, Inf)]
@parameters k [tunable = true, bounds = (0, Inf)]
eqs = [D(x) ~ (-x + k * u) / T # A first-order system with time constant T and gain k
       y ~ x]
sys = ODESystem(eqs, t, name = :tunable_first_order)
```

```@example metadata
p = tunable_parameters(sys) # extract all parameters marked as tunable
```

```@example metadata
lb, ub = getbounds(p) # operating on a vector, we get lower and upper bound vectors
```

```@example metadata
b = getbounds(sys) # Operating on the system, we get a dict
```

See also: [`ModelingToolkit.dump_variable_metadata`](@ref), [`ModelingToolkit.dump_parameters`](@ref),
[`ModelingToolkit.dump_unknowns`](@ref).

## Index

```@index
Pages = ["Variable_metadata.md"]
```

## Docstrings

```@autodocs
Modules = [ModelingToolkit]
Pages = ["variables.jl"]
Private = false
```

```@docs
ModelingToolkit.dump_variable_metadata
ModelingToolkit.dump_parameters
ModelingToolkit.dump_unknowns
```
