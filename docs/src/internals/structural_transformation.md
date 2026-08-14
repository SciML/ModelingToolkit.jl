# Structural Transformation

!!! warning "Developer API"
    These APIs support ModelingToolkit's structural-transformation machinery and
    extension packages. They are version-controlled for ModelingToolkit developers;
    end-user applications should use the documented public simplification API instead.

These functions are used for structural analysis and transformation of equation systems, including index reduction, tearing, and other algebraic manipulations used in the simplification process.

## Tearing and Algebraic Simplification

```@docs
ModelingToolkit.tearing
ModelingToolkit.tearing_substitution
ModelingToolkit.dummy_derivative
```

## Index Reduction

```@docs
ModelingToolkit.dae_index_lowering
ModelingToolkit.pantelides_reassemble
```

## Incidence Matrix Operations

```@docs
ModelingToolkit.sorted_incidence_matrix
ModelingToolkit.but_ordered_incidence
```

## Variable Ordering and Masks

```@docs
ModelingToolkit.lowest_order_variable_mask
ModelingToolkit.highest_order_variable_mask
```

## In-Place Compilation

```@docs
ModelingToolkit.mtkcompile!
```
