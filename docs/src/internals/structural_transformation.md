# Structural Transformation

!!! warning "Developer API"
    These APIs support ModelingToolkit's structural-transformation machinery and
    extension packages. They are version-controlled for ModelingToolkit developers;
    end-user applications should use the documented public simplification API instead.

These functions are used for structural analysis and transformation of equation systems, including index reduction, tearing, and other algebraic manipulations used in the simplification process.

The names below are versioned developer API. `TearingState` is mutable compiler state and
its fields are intentionally opaque. `find_solvables!` is retained as a compatibility
export; new code should call the documented transformation entry points instead.

## Tearing and Algebraic Simplification

```@docs
ModelingToolkit.tearing
ModelingToolkit.tearing_substitution
ModelingToolkit.dummy_derivative
StateSelection.find_solvables!
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
ModelingToolkitTearing.TearingState
```

## Structural Transformation Interfaces

These upstream interfaces are rendered here because the ModelingToolkit structural
transformation developer API dispatches on them.

```@docs
StateSelection.TransformationState
StateSelection.TearingAlgorithm
StateSelection.TearingResult
StateSelection.find_eq_solvables!
StateSelection.bareiss.bareiss!
StateSelection.CLIL.SparseMatrixCLIL
ModelingToolkitTearing.InlineLinearSystem
ModelingToolkitTearing.inline_linear_systems
```
