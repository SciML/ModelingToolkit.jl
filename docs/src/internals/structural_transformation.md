# Structural Transformation

!!! warning "Internal API"
    The functions documented on this page are internal implementation details of ModelingToolkit. They are not part of the public API and may change or be removed without notice in non-breaking releases. This documentation is provided to help contributors understand the codebase.

These functions are used for structural analysis and transformation of equation systems, including index reduction, tearing, and other algebraic manipulations used in the simplification process.

## Tearing and Algebraic Simplification

```@docs
tearing
tearing_reassemble
tearing_substitution
torn_system_jacobian_sparsity
find_solvables!
linear_subsys_adjmat!
```

## Index Reduction

```@docs
dae_index_lowering
pantelides!
pantelides_reassemble
dummy_derivative
```

## Consistency Checking

```@docs
check_consistency
```

## Incidence Matrix Operations

```@docs
sorted_incidence_matrix
but_ordered_incidence
```

## Variable Ordering and Masks

```@docs
lowest_order_variable_mask
highest_order_variable_mask
computed_highest_diff_variables
```

## Shift Operations

These functions handle shift operations in discrete-time systems.

```@docs
shift2term
lower_shift_varname
simplify_shifts
distribute_shift
```

## System Structure Types and Functions

```@docs
SystemStructure
TearingState
TransformationState
isdiffvar
isdervar
isalgvar
isdiffeq
algeqs
is_only_discrete
dervars_range
diffvars_range
algvars_range
get_fullvars
system_subset
```

## Graph Types

```@docs
Matching
InducedCondensationGraph
MatchedCondensationGraph
Unassigned
unassigned
```