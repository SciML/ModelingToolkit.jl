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
