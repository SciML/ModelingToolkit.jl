# API


## Intermediate representation

### Types

```@docs
Expression
Variable
ModelingToolkit.Constant
Operation
Equation
Differential
```

### Functions

```@docs
Base.get(c::ModelingToolkit.Constant)
Base.:~(::Expression, ::Expression)
expand_derivatives
ModelingToolkit.derivative
```

### Macros

```@docs
@parameters
@variables
@derivatives
@register
```

## Systems

### Types

```@docs
ModelingToolkit.AbstractSystem
ODESystem
NonlinearSystem
```

### Functions

```@docs
independent_variables
dependent_variables
parameters
calculate_jacobian
generate_jacobian
generate_function
DiffEqBase.ODEFunction(sys::ODESystem, dvs, ps; version::FunctionVersion = ArrayFunction)
```
