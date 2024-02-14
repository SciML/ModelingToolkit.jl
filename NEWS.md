# ModelingToolkit v9 Release Notes

### Upgrade guide

- The function `states` is renamed to `unknowns`. In a similar vein:
  - `unknown_states` is now `solved_unknowns`.
  - `get_states` is `get_unknowns`.
  - `get_unknown_states` is now `get_solved_unknowns`.
- The default backend for using units in models is now `DynamicQuantities.jl` instead of
  `Unitful.jl`.
- ModelingToolkit.jl now exports common definitions of `t` (time independent variable)
  and `D` (the first derivative with respect to `t`). Any models made using ModelingToolkit.jl
  should leverage these common definitions. There are three variants:
  - `t` and `D` use DynamicQuantities.jl units. This is the default for standard library
    components.
  - `t_unitful` and `D_unitful` use Unitful.jl units.
  - `t_nounits` and `D_nounits` are unitless.
- `ODAEProblem` is deprecated in favor of `ODEProblem`.
- Specifying the independent variable for an `ODESystem` is now mandatory. The `ODESystem(eqs)`
  constructor is removed.
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
  - To restore the old behavior, use `ModelingToolkit.@set sys.index_cache = nothing` before creating
    a problem, and after calling `structural_simplify`.
- Discrete-time system using `Difference` are unsupported. Instead, use the new `Clock`-based syntax.
