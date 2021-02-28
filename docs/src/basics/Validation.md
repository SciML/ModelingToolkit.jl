# Model Validation and Units

ModelingToolkit.jl provides extensive functionality for model validation
and unit checking. This is done by providing metadata to the variable
types and then running the validation functions which identify malformed
systems and non-physical equations.

## Consistency Checking

```@docs
check_consistency
```

## Unit and Type Validation

```@docs
ModelingToolkit.validate
```
