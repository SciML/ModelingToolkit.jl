# [Symbolic variables and variable metadata](@id symbolic_metadata)

ModelingToolkit uses [Symbolics.jl](https://docs.sciml.ai/Symbolics/stable/) for the symbolic
manipulation infrastructure. In fact, the `@variables` macro is defined in Symbolics.jl. In
addition to `@variables`, ModelingToolkit defines `@parameters`, `@independent_variables`,
`@constants` and `@brownians`. These macros function identically to `@variables` but allow
ModelingToolkit to attach additional metadata.

```@docs
Symbolics.@variables
@independent_variables
@parameters
@constants
@brownians
```

Symbolic variables can have metadata attached to them. The defaults and guesses assigned
at variable construction time are examples of this metadata. ModelingToolkit also defines
additional types of metadata.

## Variable defaults

Variables can be assigned default values to avoid having to specify defaults to the
[`System`](@ref) constructor.

```@docs
ModelingToolkit.hasdefault
ModelingToolkit.getdefault
ModelingToolkit.setdefault
```

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
@named sys = System([u ~ p], t)
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

```@docs
hasdescription
getdescription
ModelingToolkit.VariableDescription
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

```@docs
hasconnect
getconnect
ModelingToolkit.VariableConnectType
```

```@docs; canonical=false
Flow
Stream
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

```@docs
isinput
isoutput
ModelingToolkit.setinput
ModelingToolkit.setoutput
ModelingToolkit.VariableInput
ModelingToolkit.VariableOutput
```

## Bounds

Bounds are useful when parameters are to be optimized, or to express intervals of uncertainty.

```@repl metadata
@variables u [bounds = (-1, 1)];
hasbounds(u)
getbounds(u)
```

Bounds can also be specified for array variables. A scalar array bound is applied to each
element of the array. A bound may also be specified as an array, in which case the size of
the array must match the size of the symbolic variable.

```@repl metadata
@variables x[1:2, 1:2] [bounds = (-1, 1)];
hasbounds(x)
getbounds(x)
getbounds(x[1, 1])
getbounds(x[1:2, 1])
@variables x[1:2] [bounds = (-Inf, [1.0, Inf])];
hasbounds(x)
getbounds(x)
getbounds(x[2])
hasbounds(x[2])
```

```@docs
hasbounds
getbounds
ModelingToolkit.VariableBounds
```

## Guess

Specify an initial guess for variables of a `System`. This is used when building the
[`InitializationProblem`](@ref).

```@repl metadata
@variables u [guess = 1];
hasguess(u)
getguess(u)
```

```@docs
hasguess
getguess
```

When a system is constructed, the guesses of the involved variables are stored in a `Dict`
in the system. After this point, the guess metadata of the variable is irrelevant.

```@docs; canonical=false
guesses
```

## Mark input as a disturbance

Indicate that an input is not available for control, i.e., it's a disturbance input.

```@example metadata
@variables u [input = true, disturbance = true]
isdisturbance(u)
```

```@docs
isdisturbance
```

## Mark parameter as tunable

Indicate that a parameter can be automatically tuned by parameter optimization or automatic control tuning apps.

```@example metadata
@parameters Kp [tunable = true]
istunable(Kp)
```

```@docs
istunable
ModelingToolkit.isconstant
```

!!! note
    
    [`@constants`](@ref) is a convenient way to create `@parameters` with `tunable = false`
    metadata

## Probability distributions

A probability distribution may be associated with a parameter to indicate either
uncertainty about its value, or as a prior distribution for Bayesian optimization.

```@repl metadata
using Distributions;
d = Normal(10, 1);
@parameters m [dist = d];
hasdist(m)
getdist(m)
```

```@docs
hasdist
getdist
```

## Irreducible

A variable can be marked `irreducible` to prevent it from being moved to an
`observed` state. This forces the variable to be computed during solving so that
it can be accessed in [callbacks](@ref events)

```@example metadata
@variables important_value [irreducible = true]
isirreducible(important_value)
```

```@docs
isirreducible
ModelingToolkit.VariableIrreducible
```

## State Priority

When a model is structurally simplified, the algorithm will try to ensure that the variables with higher state priority become states of the system. A variable's state priority is a number set using the `state_priority` metadata.

```@example metadata
@variables important_dof [state_priority = 10] unimportant_dof [state_priority = -2]
state_priority(important_dof)
```

```@docs
state_priority
ModelingToolkit.VariableStatePriority
```

## Units

Units for variables can be designated using symbolic metadata. For more information, please see the [model validation and units](@ref units) section of the docs. Note that `getunit` is not equivalent to `get_unit` - the former is a metadata getter for individual variables (and is provided so the same interface function for `unit` exists like other metadata), while the latter is used to handle more general symbolic expressions.

```@repl metadata
using DynamicQuantities;
@variables speed [unit = u"m/s"];
hasunit(speed)
getunit(speed)
```

```@docs
hasunit
getunit
ModelingToolkit.VariableUnit
```

## Variable type

This metadata is used by the [`System`](@ref) constructor for automatically identifying the different types of variables in a system.

```@docs
ModelingToolkit.VariableType
ModelingToolkit.MTKVariableTypeCtx
ModelingToolkit.isparameter
```

## Miscellaneous metadata

User-defined metadata can be added using the `misc` metadata. This can be queried
using the `hasmisc` and `getmisc` functions.

```@repl metadata
@variables u [misc = :conserved_parameter] y [misc = [2, 4, 6]];
hasmisc(u)
getmisc(y)
```

```@docs
hasmisc
getmisc
ModelingToolkit.VariableMisc
```

## Dumping metadata

ModelingToolkit allows dumping the metadata of a variable as a `NamedTuple`.

```@docs
ModelingToolkit.dump_variable_metadata
```

## Additional functions

For systems that contain parameters with metadata like described above, have some additional functions defined for convenience.
In the example below, we define a system with tunable parameters and extract bounds vectors

```@example metadata
@variables x(t)=0 u(t)=0 [input=true] y(t)=0 [output=true]
@parameters T [tunable = true, bounds = (0, Inf)]
@parameters k [tunable = true, bounds = (0, Inf)]
eqs = [D(x) ~ (-x + k * u) / T # A first-order system with time constant T and gain k
       y ~ x]
sys = System(eqs, t, name = :tunable_first_order)
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

See also:

```@docs; canonical=false
tunable_parameters
ModelingToolkit.dump_unknowns
ModelingToolkit.dump_parameters
```

## Symbolic operators

ModelingToolkit makes heavy use of "operators". These are custom functions that are applied
to symbolic variables. The most common operator is the `Differential` operator, defined in
Symbolics.jl.

```@docs
Symbolics.Differential
```

ModelingToolkit also defines a plethora of custom operators.

```@docs
Pre
Initial
Shift
EvalAt
```

While not an operator, `ShiftIndex` is commonly used to use `Shift` operators in a more
convenient way when writing discrete systems.

```@docs
ShiftIndex
```

### Sampled time operators

The following operators are used in hybrid ODE systems, where part of the dynamics of the
system happen at discrete intervals on a clock. While ModelingToolkit cannot yet simulate
such systems, it has the capability to represent them.

!!! warn
    
    These operators are considered experimental API.

```@docs
Sample
Hold
SampleTime
sampletime
```
