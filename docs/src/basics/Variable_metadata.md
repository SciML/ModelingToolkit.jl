# Symbolic metadata
It is possible to add metadata to symbolic variables. The following
information can be added (note, it's possible to extend this to user-defined metadata as well)

## Input or output
Designate a variable as either an input or an output using the following
```@example metadata
using ModelingToolkit
@variables u [input=true]
isinput(u)
```
```@example metadata
@variables y [output=true]
isoutput(y)
```

## Bounds
Bounds are useful when parameters are to be optimized, or to express intervals of uncertainty.

```@example metadata
@variables u [bounds=(-1,1)]
hasbounds(u)
```
```@example metadata
getbounds(u)
```

## Mark input as a disturbance 
Indicate that an input is not available for control, i.e., it's a disturbance input.

```@example metadata
@variables u [input=true, disturbance=true]
isdisturbance(u)
```

## Mark parameter as tunable
Indicate that a parameter can be automatically tuned by automatic control tuning apps.

```@example metadata
@parameters Kp [tunable=true]
istunable(Kp)
```

## Probability distributions
A probability distribution may be associated with a parameter to indicate either
uncertainty about it's value, or as a prior distribution for Bayesian optimization.

```julia
using Distributions
d = Normal(10, 1)
@parameters m [dist=d]
hasdist(m)
```
```julia
getdist(m)
```

## Additional functions
For systems that contain parameters with metadata like described above have some additional functions defined for convenience.
In the example below, we define a system with tunable parameters and extract bounds vectors

```@example metadata
@parameters t
Dₜ = Differential(t)
@variables x(t)=0 u(t)=0 [input=true] y(t)=0 [output=true]
@parameters T [tunable = true, bounds = (0, Inf)]
@parameters k [tunable = true, bounds = (0, Inf)]
eqs = [
    Dₜ(x) ~ (-x + k*u) / T # A first-order system with time constant T and gain k
    y ~ x
]
sys = ODESystem(eqs, t, name=:tunable_first_order)
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
