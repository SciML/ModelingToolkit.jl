# ModelingToolkit v11 Release Notes

## Symbolics@7 and SymbolicUtils@4 compatibility

SymbolicUtils version 4 involved a major overhaul of the core symbolic infrastructure, which
propagated to Symbolics as Symbolics version 7. ModelingToolkit has now updated to these versions.
This includes significant type-stability improvements, enabling precompilation of large parts
of the symbolic infrastructure and faster TTFX. It is highly recommended to read the
[Release Notes for SymbolicUtils@4](https://github.com/JuliaSymbolics/SymbolicUtils.jl/releases/tag/v4.0.0)
and the [doc page](https://docs.sciml.ai/SymbolicUtils/dev/manual/variants/) describing the new
variant structure before these release notes.

As part of these changes, ModelingToolkit has changed how some data is represented to allow
precompilation. Notably, `variable => value` mappings (such as guesses) are stored as an
`AbstractDict{SymbolicT, SymbolicT}`. Here, `SymbolicT` is a type that comes from Symbolics.jl,
and is the type for all unwrapped symbolic values. This means that any non-symbolic values
are stored as `SymbolicUtils.Const` variants. Mutation such as `guesses(sys)[x] = 1.0` is still
possible, and values are automatically converted. However, obtaining the value back requires
usage of `SymbolicUtils.unwrap_const` or `Symbolics.value`.

Following is a before/after comparison of the TTFX for the most common operations in ModelingToolkit.jl.
Further improvements are ongoing. Note that the timings do depend on many factors such as the exact system
used, types passed to constructor functions, other packages currently loaded in the session, presence of
array variables/equations, whether index reduction is required, and the behavior of various passes in
`mtkcompile`. However, the numbers are good representations of the kinds of performance improvements
that are possible due to the new infrastructure. There will continue to be improvements as this gets
more extensive testing and we are better able to identify bottlenecks in compilation.

### `System` constructor

The time to call `System`, not including the time taken for `@variables` or building the equations.

Before:

```
  0.243758 seconds (563.80 k allocations: 30.613 MiB, 99.48% compilation time: 3% of which was recompilation)
elapsed time (ns):  2.43757958e8
gc time (ns):       0
bytes allocated:    32099616
pool allocs:        563137
non-pool GC allocs: 16
malloc() calls:     651
free() calls:       0
minor collections:  0
full collections:   0
```

After:

```
  0.000670 seconds (217 allocations: 10.641 KiB)
elapsed time (ns):  669875.0
gc time (ns):       0
bytes allocated:    10896
pool allocs:        217
non-pool GC allocs: 0
minor collections:  0
full collections:   0
```

### `complete`

Before:

```
  1.795140 seconds (9.76 M allocations: 506.143 MiB, 2.67% gc time, 99.75% compilation time: 71% of which was recompilation)
elapsed time (ns):  1.795140083e9
gc time (ns):       47998414
bytes allocated:    530729216
pool allocs:        9747214
non-pool GC allocs: 111
malloc() calls:     10566
free() calls:       8069
minor collections:  5
full collections:   1
```

After:

```
  0.001191 seconds (1.08 k allocations: 2.554 MiB)
elapsed time (ns):  1.190625e6
gc time (ns):       0
bytes allocated:    2678088
pool allocs:        1077
non-pool GC allocs: 0
malloc() calls:     3
free() calls:       0
minor collections:  0
full collections:   0
```

### `TearingState` constructor

`TearingState` is an intermediary step in `mtkcompile`. It is significant enough for the impact
to be worth measuring separately.

Before:

```
  0.374312 seconds (527.01 k allocations: 32.318 MiB, 24.13% gc time, 99.60% compilation time: 85% of which was recompilation)
elapsed time (ns):  3.74312e8
gc time (ns):       90318708
bytes allocated:    33888248
pool allocs:        526440
non-pool GC allocs: 11
malloc() calls:     555
free() calls:       2923
minor collections:  1
full collections:   0
```

After:

```
  0.002062 seconds (1.07 k allocations: 8.546 MiB, 50.24% compilation time)
elapsed time (ns):  2.0618339999999998e6
gc time (ns):       0
bytes allocated:    8961560
pool allocs:        1064
non-pool GC allocs: 0
malloc() calls:     6
free() calls:       0
minor collections:  0
full collections:   0
```

### `mtkcompile`

This measures the time taken by the first call to `mtkcompile`. This is run after the `TearingState`
benchmark, and hence the compile time from that aspect of the process is not included (runtime is
included).

Before:

```
  1.772756 seconds (3.81 M allocations: 206.068 MiB, 0.63% gc time, 99.71% compilation time: 71% of which was recompilation)
elapsed time (ns):  1.772755875e9
gc time (ns):       11162292
bytes allocated:    216077752
pool allocs:        3808615
non-pool GC allocs: 61
malloc() calls:     4877
free() calls:       4844
minor collections:  2
full collections:   0
```

After:

```
  0.018629 seconds (20.74 k allocations: 932.062 KiB, 89.89% compilation time)
elapsed time (ns):  1.8628542e7
gc time (ns):       0
bytes allocated:    954432
pool allocs:        20727
non-pool GC allocs: 0
malloc() calls:     13
free() calls:       0
minor collections:  0
full collections:   0
```

## Semantic separation of discretes

ModelingToolkit has long overloaded the meaning of `@parameters` to the point that it means
"anything that isn't `@variables`." This isn't a very intuitive or clear definition. This is
now improved with the introduction of `@discretes`. Any quantities that vary on a different
time-scale than those in `@variables` are now `@discretes`. `@parameters` can only be used to
create "time-independent parameters". For clarity, the following continues to work:

```julia
@parameters p q[1:3] f(::Real, ::Real)
```

However, this is now disallowed:

```julia
@parameters value(t)
```

Instead, it must be declared as:

```julia
@discretes value(t)
```

And can be passed along with the `@variables`. Essentially, for time-varying systems
the constructor syntax is

```julia
System(equations, independent_variable, time_varying_variables, constant_values, [brownians])
```

In the subsequent release notes and in documentation, "variables" refers to either `@variables`
or `@discretes` unless explicitly mentioned otherwise.

An important note is that while this is a difference in declaration, the semantics are defined
by their usage in the system. More concretely, a variable declared via `@discretes` is only
actually considered discrete if it is part of the variables updated in a callback in the system.
Otherwise, it is treated identically to a variable declared via `@variables`.

## Changes to `defaults` and initialization semantics

The concept of `defaults` is a relic of earlier ModelingToolkit versions, from when initialization
did not exist and they served as convenient initial conditions. The package has evolved greatly since then
and `defaults` have taken on many different meanings in different contexts. This makes their usage
complicated and unintuitive.

`defaults` have now been removed. They are replaced by two new concepts, with simple and well-defined
semantics. Firstly, `initial_conditions` is a variable-value mapping aimed solely at being a convenient
way to provide initial conditions to `SciMLProblem`s constructed from the system. Specifying them is
identical to providing initial values to the `ODEProblem` constructor. Secondly, `bindings` is an
immutable variable-value mapping representing strong constraints between variables/parameters.
A binding for a variable is a function of other variables/parameters that is enforced during initialization.
A binding for a parameter is a function of other parameters that exclusively defines the value of that
parameter. Bound variables or parameters cannot be given initial conditions, either through the
`initial_conditions` keyword or by passing them to the problem constructor. In effect, bindings
serve to mark specific variables as aliases of others during initialization, and parameters as aliases
of other parameters. This supersedes the previous concept of parameter bindings, and explicit parameter
equations passed along with the equations of the model. Since bound parameters are computed as functions
of other parameters, they are treated akin to observed variables. They are not stored in the parameter
object, and instead are computed on the fly as required.

Sometimes, it is useful to enforce a relation between parameters while allowing them to be given initial
values. For example, one might relate the radius `r` and area `A` of a pipe as `A ~ pi * r * r`. Users of
the model should be able to provide a value for either `r` or `A`, and the other should be calculated
automatically. This is done by providing the relation `A ~ pi * r * r` to the `initialization_eqs`
keyword of the model and binding both `A` and `r` to `missing`. Similar to v10, the equation represents
a constraint to be enforced. The bindings act similar to the `missing` defaults in v9 and v10, indicating
that the parameters are to be solved for. They are part of bindings since a parameter to be solved for
cannot be an alias for a different value. As such, the choice of parameters that can be solved for is
an immutable property of the system. Note that making a parameter solvable no longer requires specifying a
guess. If a guess is required to solve the initialization, ModelingToolkit will error with an informative
message during problem construction. Note that since parameters can only be bound to other parameters,
a parameter `x0` can be bound to the initial value of a variable `x` using the binding `x0 = Initial(x)`.

The formulation of the initialization system can now be summarized succinctly. The system solves for:

- Unknowns of the system.
- Observables (observed variables) of the system.
- All unknowns for which derivatives are known (differential variables, and ones for which derivative
  information is available due to the index reduction process).
- Discrete variables (created via `@discretes`).
- Parameters with a binding of `missing`.

It is composed of:

- Algebraic equations.
- Observed equations.
- Initialization equations.
- The `initial_conditions` of the system.
- Initial conditions passed to the problem constructor. These override values in `initial_conditions`
  for the same variable.

Additionally, `Initial` parameters exist for the following variables:

- Unknowns
- Observables
- First derivatives of all unknowns and observables
- Discrete variables
- Parameters with a binding of `missing`

"Defaults" specified via variable metadata are now translated into either `initial_conditions` or
`bindings` depending on the value. If the value is a constant, it is part of `initial_conditions`.
If it is an expression involving other variables/parameters, it is part of `bindings`. For example,
the following are `initial_conditions`:

```julia
@variables x(t) = 1 y(t)[1:3] = zeros(3)
@parameters f(::Real) = sin
```

The following are bindings:

```julia
@variables z(t) = x w(t)[1:2] = [1.5, z]
@parameters p[1:3] = f(3)
```

Notably, arrays are considered atomic. This means that if even one element of an array default is
symbolic, the entire array variable is considered bound. Partial bindings can be constructed by
destructuring the array:

```julia
@parameters par[1:3] = [par1, par2, par3]
```

Where `par1`, `par2` and `par3` can independently have initial conditions or bindings. In a
similar vein, `guesses`, `initial_conditions` and `bindings` are all stored in special
`AbstractDict` types that disallow scalarized keys. For example, `par[1]` cannot be a key
of these dictionaries. `par` is allowed as a key. Initial values can still be given to the
problem constructor in scalarized form.

As mentioned previously, bindings cannot be mutated. To change the bindings of a system,
the following pattern can be employed:

```julia
binds = bindings(sys)
# the `ReadOnlyDict` wrapper uses `Base.parent` to get the underlying mutable container
new_binds = parent(copy(binds))

# mutate `new_binds`...

using Setfield: @set!

@set! sys.bindings = new_binds
sys = complete(sys) # Important!
```

Mutation of bindings without copying them is undefined behavior and can lead to unpredictable bugs.

## Array variables as inputs

Previously, ModelingToolkit allowed part of an array variable to be an input. For example, the following
used to be valid:

```julia
@variables u(t)[1:2] [input = true]
@named sys = # Some system involving `u`

sys = mtkcompile(sys; inputs = [u[1]])
```

This is now disallowed. `mtkcompile` will throw an informative error if part of an array is passed as an
input.

## Deprecation of `@mtkmodel`

The `@mtkmodel` originated as a convenient DSL for creating models. However, it has not received the same
level of support as other features due to the complexity of the parsing. It is also a major source of bugs,
and is thus a tripping hazard for new and old users alike. The macro is now deprecated. It is moved to a new
package, SciCompDSL.jl. Enough updates have been made to allow it to create systems in v11, but it will not
receive more maintenance from the core developers. It is, however, still open to community contribution. For
more details, please refer to the discussion in [this Discourse thread](https://discourse.julialang.org/t/using-mtk-when-i-import-modelingtoolkit/133681/12).

## Splitting into `ModelingToolkitBase` and relicensing of parts of ModelingToolkit

The advanced structural simplification algorithms in ModelingToolkit, such as index reduction, structural
singularity removal and tearing, are now moved to the [StateSelection.jl](https://github.com/JuliaComputing/StateSelection.jl/)
and ModelingToolkitTearing.jl (`lib/ModelingToolkitTearing` in the same repo) packages. These packages are
AGPL licensed. ModelingToolkitBase.jl contains the `System` representation, callbacks, all of the code-generation
targets (problem constructors), and initialization infrastructure. Items that depend on structural simplification
and/or require highly specialized code generation (`SCCNonlinearProblem` prominent among them) have been moved
to ModelingToolkit.jl, which depends on the aforementioned AGPL packages. ModelingToolkitBase still contains
a simple version of `mtkcompile` suitable for most use cases. It does not perform index reduction and requires
that all differential equations are explicit in the derivative. However, it does have a simpler tearing algorithm
that is capable of identifying observed equations in many scenarios. In fact, it is also able to reduce many
systems created using modular components and the `connect` infrastructure. Contributions to improve this
are also welcome.

For more information on the split and surrounding changes, please follow the discussion in
[this Discourse thread](https://discourse.julialang.org/t/modelingtoolkit-v11-library-split-and-licensing-community-feedback-requested/134396).

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
