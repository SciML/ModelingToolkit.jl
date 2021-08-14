# [Model Validation and Units](@id units)

ModelingToolkit.jl provides extensive functionality for model validation and unit checking. This is done by providing metadata to the variable types and then running the validation functions which identify malformed systems and non-physical equations. This approach provides high performance and compatibility with numerical solvers.

## Assigning Units 

Units may assigned with the following syntax. 

```julia
using ModelingToolkit, Unitful
@variables t [unit = u"s"] x(t) [unit = u"m"] g(t) w(t) [unit = "Hz"]
#Or,
@variables(t, [unit = u"s"], x(t), [unit = u"m"], g(t), w(t), [unit = "Hz"])
#Or,
@variables(begin
t, [unit = u"s"],
x(t), [unit = u"m"],
g(t),
w(t), [unit = "Hz"]
end)
```

Do not use `quantities` such as `1u"s"` or `1/u"s"` or `u"1/s"` as these will result in errors; instead use `u"s"` or `u"s^1"`. 

## Unit Validation & Inspection

Unit validation of equations happens automatically when creating a system. However, for debugging purposes one may wish to validate the equations directly using `validate`.

```@docs
ModelingToolkit.validate
```

Inside, `validate` uses `get_unit`, which may be directly applied to any term. Note that `validate` will not throw an error in the event of incompatible units, but `get_unit` will. If you would rather receive a warning instead of an error, use `safe_get_unit` which will yield `nothing` in the event of an error. Unit agreement is tested with `ModelingToolkit.equivalent(u1,u2)`. 


```@docs
ModelingToolkit.get_unit
```

Example usage below. Note that `ModelingToolkit` does not force unit conversions to preferred units in the event of nonstandard combinations -- it merely checks that the equations are consistent. 

```julia
using ModelingToolkit, Unitful
@parameters τ [unit = u"ms"]
@variables t [unit = u"ms"] E(t) [unit = u"kJ"] P(t) [unit = u"MW"]
D = Differential(t)
eqs = eqs = [D(E) ~ P - E/τ,
                0 ~ P       ]
ModelingToolkit.validate(eqs) #Returns true
ModelingToolkit.validate(eqs[1]) #Returns true
ModelingToolkit.get_unit(eqs[1].rhs) #Returns u"kJ ms^-1"
```

An example of an inconsistent system: at present, `ModelingToolkit` requires that the units of all terms in an equation or sum to be equal-valued (`ModelingToolkit.equivalent(u1,u2)`), rather that simply dimensionally consistent. In the future, the validation stage may be upgraded to support the insertion of conversion factors into the equations. 

```julia
using ModelingToolkit, Unitful
@parameters τ [unit = u"ms"]
@variables t [unit = u"ms"] E(t) [unit = u"J"] P(t) [unit = u"MW"]
D = Differential(t)
eqs = eqs = [D(E) ~ P - E/τ,
                0 ~ P       ]
ModelingToolkit.validate(eqs) #Returns false while displaying a warning message
```

## `Unitful` Literals & User-Defined Functions

In order for a function to work correctly during both validation & execution, the function must be unit-agnostic. That is, no unitful literals may be used. Any unitful quantity must either be a `parameter` or `variable`. For example, these equations will not validate successfully. 

```julia
using ModelingToolkit, Unitful
@variables t [unit = u"ms"] E(t) [unit = u"J"] P(t) [unit = u"MW"]
D = Differential(t)
eqs = [D(E) ~ P - E/1u"ms"   ]
ModelingToolkit.validate(eqs) #Returns false while displaying a warning message

myfunc(E) = E/1u"ms"
eqs = [D(E) ~ P - myfunc(E) ]
ModelingToolkit.validate(eqs) #Returns false while displaying a warning message
```

Instead, they should be parameterized:

```julia
using ModelingToolkit, Unitful
@parameters τ [unit = u"ms"]
@variables t [unit = u"ms"] E(t) [unit = u"kJ"] P(t) [unit = u"MW"]
D = Differential(t)
eqs = [D(E) ~ P - E/τ]
ModelingToolkit.validate(eqs) #Returns true

myfunc(E,τ) = E/τ 
eqs = [D(E) ~ P - myfunc(E,τ)]
ModelingToolkit.validate(eqs) #Returns true
```

It is recommended *not* to circumvent unit validation by specializing user-defined functions on `Unitful` arguments vs. `Numbers`. This both fails to take advantage of `validate` for ensuring correctness, and may cause in errors in the
future when `ModelingToolkit` is extended to support eliminating `Unitful` literals from functions.

## Other Restrictions

`Unitful` provides non-scalar units such as `dBm`, `°C`, etc. At this time, `ModelingToolkit` only supports scalar quantities. Additionally, angular degrees (`°`) are not supported because trigonometric functions will treat plain numerical values as radians, which would lead systems validated using degrees to behave erroneously when being solved. 

## Troubleshooting & Gotchas

If a system fails to validate due to unit issues, at least one warning message will appear, including a line number as well as the unit types and expressions that were in conflict. Some system constructors re-order equations before the unit checking can be done, in which case the equation numbers may be inaccurate. The printed expression that the problem resides in is always correctly shown.

Symbolic exponents for unitful variables *are* supported (ex: `P^γ` in thermodynamics). However, this means that `ModelingToolkit` cannot reduce such expressions to `Unitful.Unitlike` subtypes at validation time because the exponent value is not available. In this case `ModelingToolkit.get_unit` is type-unstable, yielding a symbolic result, which can still be checked for symbolic equality with `ModelingToolkit.equivalent`. 

## Parameter & Initial Condition Values

Parameter and initial condition values are supplied to problem constructors as plain numbers, with the understanding that they have been converted to the appropriate units. This is done for simplicity of interfacing with optimization solvers. Some helper function for dealing with value maps:

```julia
remove_units(p::Dict) = Dict(k => Unitful.ustrip(ModelingToolkit.get_unit(k),v) for (k,v) in p)
add_units(p::Dict) = Dict(k => v*ModelingToolkit.get_unit(k) for (k,v) in p)
```

Recommended usage:

```julia
pars = @parameters τ [unit = u"ms"]
p = Dict(τ => 1u"ms")
ODEProblem(sys,remove_units(u0),tspan,remove_units(p))
```
