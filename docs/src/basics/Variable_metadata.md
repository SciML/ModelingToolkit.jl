# [Symbolic Metadata](@id symbolic_metadata)

It is possible to add metadata to symbolic variables, the metadata will be displayed when calling help on a variable.

The following information can be added (note, it's possible to extend this to user-defined metadata as well)

## Variable descriptions

Descriptive strings can be attached to variables using the `[description = "descriptive string"]` syntax:

```@example metadata
using ModelingToolkit
@variables u [description = "This is my input"]
getdescription(u)
```

When variables with descriptions are present in systems, they will be printed when the system is shown in the terminal:

```@example metadata
@parameters t
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

@variables t, i(t) [connect = Flow]
@variables k(t) [connect = Stream]
```

## Input or output

Designate a variable as either an input or an output using the following

```@example metadata
using ModelingToolkit
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

## Additional functions

For systems that contain parameters with metadata like described above, have some additional functions defined for convenience.
In the example below, we define a system with tunable parameters and extract bounds vectors

```@example metadata
@parameters t
Dₜ = Differential(t)
@variables x(t)=0 u(t)=0 [input = true] y(t)=0 [output = true]
@parameters T [tunable = true, bounds = (0, Inf)]
@parameters k [tunable = true, bounds = (0, Inf)]
eqs = [Dₜ(x) ~ (-x + k * u) / T # A first-order system with time constant T and gain k
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
