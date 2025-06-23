# ModelingToolkit v10 Release Notes

## Callbacks

Callback semantics have changed.

  - There is a new `Pre` operator that is used to specify which values are before the callback.
    For example, the affect `A ~ A + 1` should now be written as `A ~ Pre(A) + 1`. This is
    **required** to be specified - `A ~ A + 1` will now be interpreted as an equation to be
    satisfied after the callback (and will thus error since it is unsatisfiable).

  - All parameters that are changed by a callback must be declared as discrete parameters to
    the callback constructor, using the `discrete_parameters` keyword argument.

```julia
event = SymbolicDiscreteCallback(
    [t == 1] => [p ~ Pre(p) + 1], discrete_parameters = [p])
```

## New `mtkcompile` and `@mtkcompile`

`structural_simplify` is now renamed to `mtkcompile`. `@mtkbuild` is renamed to
`@mtkcompile`. Their functionality remains the same. However, instead of a second
positional argument `structural_simplify(sys, (inputs, outputs))` the inputs and outputs
should be specified via keyword arguments as `mtkcompile(sys; inputs, outputs, disturbance_inputs)`.

## Reduce reliance on metadata in `mtkcompile`

Previously, `mtkcompile` (formerly `structural_simplify`) used to rely on the metadata of
symbolic variables to identify variables/parameters/brownians. This was regardless of
what the system expected the variable to be. Now, it respects the information in the system.

## Unified `System` type

There is now a single common `System` type for all types of models except PDEs, for which
`PDESystem` still exists. It follows the same syntax as `ODESystem` and `NonlinearSystem`
did. `System(equations, t[, vars, pars])` will construct a time-dependent system.
`System(equations[, vars, pars])` will construct a time-independent system. Refer to the
docstring for `System` for further information.

Utility constructors are defined for:

  - `NonlinearSystem(sys)` to convert a time-dependent system to a time-independent one for
    its steady state.
  - `SDESystem(sys, noise_eqs)` to add noise to a system
  - `JumpSystem(jumps, ...)` to define a system with jumps. Note that normal equations can
    also be passed to `jumps`.
  - `OptimizationSystem(cost, ...)` to define a system for optimization.

All problem constructors validate that the system matches the expected structure for
that problem.

## No more `parameter_dependencies`

The `parameter_dependencies` keyword is deprecated. All equations previously passed here
should now be provided as part of the standard equations of the system. If passing parameters
explicitly to the `System` constructor, the dependent parameters (on the left hand side of
parameter dependencies) should also be provided. These will be separated out when calling
`complete` or `mtkcompile`. Calling `parameter_dependencies` or `dependent_parameters` now
requires that the system is completed. The new `SDESystem` constructor still retains the
`parameter_dependencies` keyword argument since the number of equations has to match the
number of columns in `noise_eqs`.

ModelingToolkit now has discretion of what parameters are eliminated using the parameter
equations during `complete` or `mtkcompile`.

## New problem and constructors

Instead of `XProblem(sys, u0map, tspan, pmap)` for time-dependent problems and
`XProblem(sys, u0map, pmap)` for time-independent problems, the syntax has changed to
`XProblem(sys, op, tspan)` and `XProblem(sys, op)` respectively. `op` refers to the
operating point, and is a variable-value mapping containing both unknowns and parameters.

`XFunction` constructors also no longer accept the list of unknowns and parameters as
positional arguments.

## Removed `DelayParentScope`

The outdated `DelayParentScope` has been removed.

## Removed `XProblemExpr` and `XFunctionExpr`

The old `XProblemExpr` and `XFunctionExpr` constructors used to build an `Expr` that
constructs `XProblem` and `XFunction` respectively are now removed. This functionality
is now available by passing `expression = Val{true}` to any problem or function constructor.

## Renaming of `generate_*` and `calculate_*` methods

Several `generate_*` methods have been renamed, along with some `calculate_*` methods.
The `generate_*` methods also no longer accept a list of unknowns and/or parameters. Refer
to the documentation for more information.

## New behavior of `getproperty` and `setproperty!`

Using `getproperty` to access fields of a system has been deprecated for a long time, and
this functionality is now removed. `setproperty!` previously used to update the default
of the accessed symbolic variable. This is not supported anymore. Defaults can be updated by
mutating `ModelingToolkit.get_defaults(sys)`.

## New behavior of `@constants`

`@constants` now creates parameters with the `tunable = false` metadata by default.

## Removed `FunctionalAffect`

`FunctionalAffect` is now removed in favor of the new `ImperativeAffect`. Refer to the
documentation for more information.

## Improved system metadata

Instead of an empty field that can contain arbitrary data, the `System` type stores metadata
identically to `SymbolicUtils.BasicSymbolic`. Metadata is stored in an immutable dictionary
keyed by a user-provided `DataType` and containing arbitrary values. `System` supports the
same `SymbolicUtils.getmetadata` and `SymbolicUtils.setmetadata` API as symbolic variables.
Refer to the documentation of `System` and the aforementioned functions for more information.

## Moved `connect` and `Connector` to ModelingToolkit

Previously ModelingToolkit used the `connect` function and `Connector` type defined in
Symbolics.jl. These have now been moved to ModelingToolkit along with the experimental
state machine API. If you imported them from Symbolics.jl, it is recommended to import from
ModelingToolkit instead.

## Always wrap with `ParentScope` in `@named`

When creating a system using `@named`, any symbolic quantities passed as keyword arguments
to the subsystem are wrapped in `ParentScope`. Previously, this would only happen if the
variable wasn't already wrapped in a `ParentScope`. However, the old behavior had issues
when passing symbolic quantities down multiple levels of the hierarchy. The `@named` macro
now always performs this wrapping.

# ModelingToolkit v9 Release Notes

### Upgrade guide

  - The function `states` is renamed to `unknowns`. In a similar vein:
    
      + `unknown_states` is now `solved_unknowns`.
      + `get_states` is `get_unknowns`.
      + `get_unknown_states` is now `get_solved_unknowns`.

  - The default backend for using units in models is now `DynamicQuantities.jl` instead of
    `Unitful.jl`.
  - ModelingToolkit.jl now exports common definitions of `t` (time independent variable)
    and `D` (the first derivative with respect to `t`). Any models made using ModelingToolkit.jl
    should leverage these common definitions. There are three variants:
    
      + `t` and `D` use DynamicQuantities.jl units. This is the default for standard library
        components.
      + `t_unitful` and `D_unitful` use Unitful.jl units.
      + `t_nounits` and `D_nounits` are unitless.
  - `ODAEProblem` is deprecated in favor of `ODEProblem`.
  - Specifying the independent variable for an `ODESystem` is now mandatory. The `ODESystem(eqs)`
    constructor is removed. Use `ODESystem(eqs,t)` instead.
  - Systems must be marked as `complete` before creating `*Function`/`*FunctionExpr`/`*Problem`/
    `*ProblemExpr`. Typically this involved using `@mtkbuild` to create the system or calling
    `structural_simplify` on an existing system.
  - All systems will perform parameter splitting by default. Problems created using ModelingToolkit.jl
    systems will have a custom struct instead of a `Vector` of parameters. The internals of this
    type are undocumented and subject to change without notice or a breaking release. Parameter values
    can be queried, updated or manipulated using SciMLStructures.jl or SymbolicIndexingInterface.jl.
    This also requires that the symbolic type of a parameter match its assigned value. For example,
    `@parameters p` will always use a `Float64` value for `p`. To use `Int` instead, use
    `@parameters p::Int`. Array-valued parameters must be array symbolics; `@parameters p = [1.0, 2.0]`
    is now invalid and must be changed to `@parameters p[1:2] = [1.0, 2.0]`. The index of a parameter
    in the system is also not guaranteed to be an `Int`, and will instead be a custom undocumented type.
    Parameters that have a default value depending on other parameters are now treated as dependent
    parameters. Their value cannot be modified directly. Whenever a parameter value is changed, dependent
    parameter values are recalculated. For example, if `@parameters p1 p2 = 3p1` then `p2` can not be
    modified directly. If `p1` is changed, then `p2` will be updated accordingly. To restore the old behavior:
    
      + Pass the `split = false` keyword to `structural_simplify`. E.g. `ss = structural_simplify(sys; split = false)`.
      + Pass `split = false` to `@mtkbuild`. E.g. `@mtkbuild sys = ODESystem(...) split = false`.
  - Discrete-time system using `Difference` are unsupported. Instead, use the new `Clock`-based syntax.
  - Automatic scalarization has been removed, meaning that vector variables need to be treated with proper vector
    equations. For example, `[p[1] => 1.0, p[2] => 2.0]` is no longer allowed in default equations, use
    `[p => [1.0, 2.0]]` instead. Also, array equations like for `@variables u[1:2]` have `D(u) ~ A*u` as an
    array equation. If the scalarized version is desired, use `scalarize(u)`.
  - Parameter dependencies are now supported. They can be specified using the syntax
    `(single_parameter => expression_involving_other_parameters)` and a `Vector` of these can be passed to
    the `parameter_dependencies` keyword argument of `ODESystem`, `SDESystem` and `JumpSystem`. The dependent
    parameters are updated whenever other parameters are modified, e.g. in callbacks.
  - Support for `IfElse.jl` has been dropped. `Base.ifelse` can be used instead.
  - DAE initialization and the solving for consistent initial conditions has been changed to use a customized
    initialization solve. This change adds `guess` semantics which are clearly delinated from the behavior of
    the defaults, where `default` (and `u0`) is designed to be always satisfied and error if unsatisfiable,
    while `guess` is an initial guess to the initializer. In previous iterations, initialization with the
    default (`BrownBasicInit`) would treat the initial condition to the algebraic variables as a `guess`,
    and with `ShampineCollocationInit` would treat all initial conditions as a `guess`. To return to the
    previous behavior, use the keyword argument `initializealg` in the solve, i.e.
    `solve(prob;initializealg = BrownBasicInit())`.
