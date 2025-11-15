# [The `System` type](@id System_type)

ModelingToolkit.jl uses `System` to symbolically represent all types of numerical problems.
Users create `System`s representing the problem they want to solve and `mtkcompile` transforms
them into a format ModelingToolkit.jl can generate code for (alongside performing other
optimizations).

```@docs
System
ModelingToolkit.AbstractSystem
```

## Utility constructors

Several utility constructors also exist to easily construct alternative system formulations.

```@docs
NonlinearSystem
SDESystem
JumpSystem
OptimizationSystem
```

## Accessor functions

Several accessor functions exist to query systems for the information they contain. In general,
for every field `x` there exists a `has_x` function which checks if the system contains the
field and a `get_x` function for obtaining the value in the field. Note that fields of a system
cannot be accessed via `getproperty` - that is reserved for accessing variables, subsystems
or analysis points of the hierarchical system.

```@docs
ModelingToolkit.has_eqs
ModelingToolkit.get_eqs
equations
ModelingToolkit.equations_toplevel
full_equations
ModelingToolkit.has_noise_eqs
ModelingToolkit.get_noise_eqs
ModelingToolkit.has_jumps
ModelingToolkit.get_jumps
jumps
ModelingToolkit.has_constraints
ModelingToolkit.get_constraints
constraints
ModelingToolkit.has_costs
ModelingToolkit.get_costs
cost
ModelingToolkit.has_consolidate
ModelingToolkit.get_consolidate
ModelingToolkit.has_unknowns
ModelingToolkit.get_unknowns
unknowns
ModelingToolkit.unknowns_toplevel
ModelingToolkit.has_ps
ModelingToolkit.get_ps
parameters
ModelingToolkit.parameters_toplevel
tunable_parameters
ModelingToolkit.has_brownians
ModelingToolkit.get_brownians
brownians
ModelingToolkit.has_iv
ModelingToolkit.get_iv
ModelingToolkit.has_observed
ModelingToolkit.get_observed
observed
observables
ModelingToolkit.has_name
ModelingToolkit.get_name
nameof
ModelingToolkit.has_description
ModelingToolkit.get_description
ModelingToolkit.description
ModelingToolkit.has_defaults
ModelingToolkit.get_defaults
defaults
ModelingToolkit.has_guesses
ModelingToolkit.get_guesses
guesses
ModelingToolkit.get_systems
ModelingToolkit.has_initialization_eqs
ModelingToolkit.get_initialization_eqs
initialization_equations
ModelingToolkit.has_continuous_events
ModelingToolkit.get_continuous_events
continuous_events
ModelingToolkit.continuous_events_toplevel
ModelingToolkit.has_discrete_events
ModelingToolkit.get_discrete_events
ModelingToolkit.discrete_events_toplevel
ModelingToolkit.has_assertions
ModelingToolkit.get_assertions
ModelingToolkit.assertions
ModelingToolkit.has_metadata
ModelingToolkit.get_metadata
SymbolicUtils.getmetadata(::ModelingToolkit.AbstractSystem, ::DataType, ::Any)
SymbolicUtils.setmetadata(::ModelingToolkit.AbstractSystem, ::DataType, ::Any)
ModelingToolkit.has_is_dde
ModelingToolkit.get_is_dde
ModelingToolkit.is_dde
ModelingToolkit.has_tstops
ModelingToolkit.get_tstops
ModelingToolkit.symbolic_tstops
ModelingToolkit.has_tearing_state
ModelingToolkit.get_tearing_state
ModelingToolkit.does_namespacing
toggle_namespacing
ModelingToolkit.iscomplete
ModelingToolkit.has_preface
ModelingToolkit.get_preface
ModelingToolkit.preface
ModelingToolkit.has_parent
ModelingToolkit.get_parent
ModelingToolkit.has_initializesystem
ModelingToolkit.get_initializesystem
ModelingToolkit.is_initializesystem
```

## `getproperty` syntax

ModelingToolkit allows obtaining in a system using `getproperty`. For a system `sys` with a
subcomponent `inner` containing variable `var`, `sys.inner.var` will obtain the appropriately
namespaced version of `var`. Note that this can also be used to access subsystems (`sys.inner`)
or analysis points.

!!! note
    
    By default, top-level systems not marked as `complete` will apply their namespace. Systems
    marked as `complete` will not do this namespacing. This namespacing behavior can be toggled
    independently of whether the system is completed using [`toggle_namespacing`](@ref) and the
    current namespacing behavior can be queried via [`ModelingToolkit.does_namespacing`](@ref).

```@docs
Base.getproperty(::ModelingToolkit.AbstractSystem, ::Symbol)
```

## Functions for querying system equations

```@docs
has_diff_eqs
has_alg_eqs
get_diff_eqs
get_alg_eqs
has_diff_equations
has_alg_equations
diff_equations
alg_equations
ModelingToolkit.is_alg_equation
ModelingToolkit.is_diff_equation
```

## String parsing

ModelingToolkit can parse system variables from strings.

```@docs
ModelingToolkit.parse_variable
```

## Dumping system data

```@docs
ModelingToolkit.dump_unknowns
ModelingToolkit.dump_parameters
```

```@docs; canonical = false
ModelingToolkit.dump_variable_metadata
```

## Inputs and outputs

```@docs
ModelingToolkit.inputs
ModelingToolkit.outputs
ModelingToolkit.bound_inputs
ModelingToolkit.unbound_inputs
ModelingToolkit.bound_outputs
ModelingToolkit.unbound_outputs
ModelingToolkit.is_bound
```

## Debugging utilities

```@docs
debug_system
```

## Input validation

The following values can be passed to the `check` keyword of `System` to toggle validation
of input. Flags can be combined with bitwise `|` and `&`.

```@docs
ModelingToolkit.CheckAll
ModelingToolkit.CheckNone
ModelingToolkit.CheckComponents
ModelingToolkit.CheckUnits
```

These can also be used by custom `AbstractSystem` subtypes.

## Utility functions

These utility functions can be useful when manipulating systems, especially when building
custom `AbstractSystem` subtypes.

```@docs
ModelingToolkit.collect_scoped_vars!
ModelingToolkit.collect_var_to_name!
ModelingToolkit.collect_vars!
ModelingToolkit.eqtype_supports_collect_vars
ModelingToolkit.modified_unknowns!
```

## Namespace manipulation

ModelingToolkit namespaces variables from subsystems when using them in a parent system to
disambiguate from identically named variables in other subsystems or the parent system. The
following functions are useful for manipulating namespacing functionality.

```@docs
ModelingToolkit.renamespace
ModelingToolkit.namespace_equations
```

## Linearization and Analysis

Functions for linearization and analysis of systems.

```@docs
linearization_ap_transform
get_sensitivity_function
get_comp_sensitivity_function
get_looptransfer_function
get_sensitivity
get_comp_sensitivity
get_looptransfer
open_loop
```

## Additional Equation Classification

```@docs
alg_equations
diff_equations
has_alg_equations
has_diff_equations
```
