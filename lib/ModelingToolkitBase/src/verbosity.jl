@verbosity_specifier MTKVerbosity begin
    sub_specifiers = (:initialization_verbosity,)

    toggles = (
        # compilation
        :state_priority_tie, :underconstrained_variables,
        :if_lifting_condition_grammar, :observed_equation_cycle,
        # initialization
        :singular_initialization, :overdetermined_initialization,
        :underdetermined_initialization, :scc_initialization_unavailable,
        # problem construction
        :cyclic_dependency, :overdetermined_constraints,
        :missing_scc_schedule, :dynamic_opt_time_grid,
        :imperative_affect_specification, :parameter_not_in_system,
        # analysis
        :empty_operating_point, :initialization_analysis,
        :no_unbound_inputs, :analysis_point_causality,
        # system construction
        :independent_variable_not_parameter, :constraint_variable_without_argument,
        :unbalanced_connector, :affect_default_independent_variable,
        :substitute_skips_events,
    )

    presets = (
        None = (
            initialization_verbosity = None(),
            state_priority_tie = Silent,
            underconstrained_variables = Silent,
            if_lifting_condition_grammar = Silent,
            observed_equation_cycle = Silent,
            singular_initialization = Silent,
            overdetermined_initialization = Silent,
            underdetermined_initialization = Silent,
            scc_initialization_unavailable = Silent,
            cyclic_dependency = Silent,
            overdetermined_constraints = Silent,
            missing_scc_schedule = Silent,
            dynamic_opt_time_grid = Silent,
            empty_operating_point = Silent,
            initialization_analysis = Silent,
            no_unbound_inputs = Silent,
            analysis_point_causality = Silent,
            imperative_affect_specification = Silent,
            parameter_not_in_system = Silent,
            independent_variable_not_parameter = Silent,
            constraint_variable_without_argument = Silent,
            unbalanced_connector = Silent,
            affect_default_independent_variable = Silent,
            substitute_skips_events = Silent,
        ),
        Minimal = (
            initialization_verbosity = Minimal(),
            state_priority_tie = WarnLevel,
            underconstrained_variables = Silent,
            if_lifting_condition_grammar = WarnLevel,
            observed_equation_cycle = Silent,
            singular_initialization = WarnLevel,
            overdetermined_initialization = Silent,
            underdetermined_initialization = Silent,
            scc_initialization_unavailable = Silent,
            cyclic_dependency = Silent,
            overdetermined_constraints = Silent,
            missing_scc_schedule = WarnLevel,
            dynamic_opt_time_grid = WarnLevel,
            empty_operating_point = Silent,
            initialization_analysis = InfoLevel,
            no_unbound_inputs = Silent,
            analysis_point_causality = WarnLevel,
            imperative_affect_specification = WarnLevel,
            parameter_not_in_system = WarnLevel,
            independent_variable_not_parameter = Silent,
            constraint_variable_without_argument = Silent,
            unbalanced_connector = Silent,
            affect_default_independent_variable = Silent,
            substitute_skips_events = WarnLevel,
        ),
        Standard = (
            initialization_verbosity = Minimal(),
            state_priority_tie = WarnLevel,
            underconstrained_variables = Silent,
            if_lifting_condition_grammar = WarnLevel,
            observed_equation_cycle = Silent,
            singular_initialization = WarnLevel,
            overdetermined_initialization = WarnLevel,
            underdetermined_initialization = WarnLevel,
            scc_initialization_unavailable = WarnLevel,
            cyclic_dependency = Silent,
            overdetermined_constraints = WarnLevel,
            missing_scc_schedule = WarnLevel,
            dynamic_opt_time_grid = WarnLevel,
            empty_operating_point = WarnLevel,
            initialization_analysis = InfoLevel,
            no_unbound_inputs = WarnLevel,
            analysis_point_causality = WarnLevel,
            imperative_affect_specification = WarnLevel,
            parameter_not_in_system = WarnLevel,
            independent_variable_not_parameter = WarnLevel,
            constraint_variable_without_argument = WarnLevel,
            unbalanced_connector = WarnLevel,
            affect_default_independent_variable = WarnLevel,
            substitute_skips_events = WarnLevel,
        ),
        Detailed = (
            initialization_verbosity = Minimal(),
            state_priority_tie = WarnLevel,
            underconstrained_variables = InfoLevel,
            if_lifting_condition_grammar = WarnLevel,
            observed_equation_cycle = InfoLevel,
            singular_initialization = WarnLevel,
            overdetermined_initialization = WarnLevel,
            underdetermined_initialization = WarnLevel,
            scc_initialization_unavailable = WarnLevel,
            cyclic_dependency = WarnLevel,
            overdetermined_constraints = WarnLevel,
            missing_scc_schedule = WarnLevel,
            dynamic_opt_time_grid = WarnLevel,
            empty_operating_point = WarnLevel,
            initialization_analysis = InfoLevel,
            no_unbound_inputs = WarnLevel,
            analysis_point_causality = WarnLevel,
            imperative_affect_specification = WarnLevel,
            parameter_not_in_system = WarnLevel,
            independent_variable_not_parameter = WarnLevel,
            constraint_variable_without_argument = WarnLevel,
            unbalanced_connector = WarnLevel,
            affect_default_independent_variable = WarnLevel,
            substitute_skips_events = WarnLevel,
        ),
        All = (
            initialization_verbosity = Minimal(),
            state_priority_tie = WarnLevel,
            underconstrained_variables = InfoLevel,
            if_lifting_condition_grammar = WarnLevel,
            observed_equation_cycle = InfoLevel,
            singular_initialization = WarnLevel,
            overdetermined_initialization = WarnLevel,
            underdetermined_initialization = WarnLevel,
            scc_initialization_unavailable = WarnLevel,
            cyclic_dependency = WarnLevel,
            overdetermined_constraints = WarnLevel,
            missing_scc_schedule = WarnLevel,
            dynamic_opt_time_grid = WarnLevel,
            empty_operating_point = WarnLevel,
            initialization_analysis = InfoLevel,
            no_unbound_inputs = WarnLevel,
            analysis_point_causality = WarnLevel,
            imperative_affect_specification = WarnLevel,
            parameter_not_in_system = WarnLevel,
            independent_variable_not_parameter = WarnLevel,
            constraint_variable_without_argument = WarnLevel,
            unbalanced_connector = WarnLevel,
            affect_default_independent_variable = WarnLevel,
            substitute_skips_events = WarnLevel,
        ),
    )

    groups = (
        compilation = (
            :state_priority_tie, :underconstrained_variables,
            :if_lifting_condition_grammar, :observed_equation_cycle,
        ),
        initialization = (
            :singular_initialization, :overdetermined_initialization,
            :underdetermined_initialization, :scc_initialization_unavailable,
        ),
        problem_construction = (
            :cyclic_dependency, :overdetermined_constraints,
            :missing_scc_schedule, :dynamic_opt_time_grid,
            :imperative_affect_specification, :parameter_not_in_system,
        ),
        analysis = (
            :empty_operating_point, :initialization_analysis,
            :no_unbound_inputs, :analysis_point_causality,
        ),
        system_construction = (
            :independent_variable_not_parameter, :constraint_variable_without_argument,
            :unbalanced_connector, :affect_default_independent_variable,
            :substitute_skips_events,
        ),
    )
end

@doc """
    MTKVerbosity(; preset = nothing, kwargs...)
    MTKVerbosity(preset::SciMLLogging.AbstractVerbosityPreset)

A [SciMLLogging.jl](https://github.com/SciML/SciMLLogging.jl) verbosity specifier
controlling diagnostic output from ModelingToolkit: the [`mtkcompile`](@ref) pipeline,
problem construction and initialization, and analysis utilities such as
`linearization_function` and `analyze_initialization_jacobian`.

Pass it via the `verbose` keyword. `mtkcompile` and the analysis utilities consume it
directly. `SciMLBase.*Problem` constructors share the `verbose` keyword with the solver
(problem keyword arguments are forwarded to `solve`), routed by type:

- an `MTKVerbosity` is consumed by ModelingToolkit and **not** placed in `prob.kwargs`;
- a `SciMLLogging` preset (`None()`, `Minimal()`, `Standard()`, `Detailed()`, `All()`)
  sets ModelingToolkit's verbosity **and** is forwarded to the solver, which accepts
  presets for its own verbosity specifier;
- a `Bool` maps to `Standard()`/`None()` for ModelingToolkit and is forwarded to the
  solver unchanged;
- any other value (e.g. a solver-specific specifier such as `DEVerbosity`) leaves
  ModelingToolkit at its default and is forwarded unchanged.

Each toggle is a `SciMLLogging.MessageLevel` (`Silent`, `DebugLevel`, `InfoLevel`,
`WarnLevel`, `ErrorLevel`). Note that raising a toggle to `ErrorLevel` makes the
corresponding message throw an error.

# Toggles

Compilation (group `compilation`):

- `state_priority_tie`: multiple variables in an alias group are tied for the highest
  `state_priority` during alias elimination. `WarnLevel` by default.
- `underconstrained_variables`: report variables found to be underconstrained during
  alias elimination. `Silent` by default; `InfoLevel` at `Detailed`/`All`.
- `if_lifting_condition_grammar`: the `IfLifting` pass encountered a condition outside
  the limited conditional grammar it supports and skipped it. `WarnLevel` by default.
- `observed_equation_cycle`: print the smallest cycle among observed equations when
  topological sorting fails (the failure itself is always an error). `Silent` by
  default; `InfoLevel` at `Detailed`/`All`.

Initialization (group `initialization`):

- `singular_initialization`: the initialization system is structurally singular.
  `WarnLevel` by default.
- `overdetermined_initialization` / `underdetermined_initialization`: the initialization
  system is over/underdetermined and will be solved in a least-squares sense.
  `WarnLevel` by default.
- `scc_initialization_unavailable`: `SCCNonlinearProblem`-based initialization requires
  `split = true`. `WarnLevel` by default.

Problem construction (group `problem_construction`):

- `cyclic_dependency`: report cycles among the initial conditions of unknowns or
  parameters. `Silent` by default; `WarnLevel` at `Detailed`/`All`.
- `overdetermined_constraints`: a `BVProblem` or dynamic optimization problem has more
  constraints and initial conditions than states. `WarnLevel` by default.
- `missing_scc_schedule`: a simplified system unexpectedly lacks a tearing schedule
  when constructing an `SCCNonlinearProblem`. `WarnLevel` by default.
- `dynamic_opt_time_grid`: a `dt` or `steps` argument to a dynamic optimization problem
  is ignored for the given time span. `WarnLevel` by default.
- `imperative_affect_specification`: an `ImperativeAffect` of a callback has a duplicate
  name, an expression that cannot be modified, or a symbol that is both observed and
  modified. `WarnLevel` by default.
- `parameter_not_in_system`: `remake` (or another `remake_buffer` call) got a symbolic
  variable that is not a parameter of the system, and ignored it. `WarnLevel` by default.

Analysis (group `analysis`):

- `empty_operating_point`: an empty operating point was passed to
  `linearization_function`. `WarnLevel` by default.
- `initialization_analysis`: the report of `analyze_initialization_jacobian`.
  `InfoLevel` by default.
- `no_unbound_inputs`: `generate_control_function` found no unbound inputs.
  `WarnLevel` by default.
- `analysis_point_causality`: an analysis-point `connect` looks like it has reversed
  causality (input where an output is expected, or vice versa). `WarnLevel` by default.

System construction (group `system_construction`):

- `independent_variable_not_parameter`: the independent variable of a `System` is not
  defined with `@independent_variables`. `WarnLevel` by default.
- `constraint_variable_without_argument`: a constraint or cost uses an unknown with no
  argument, so it applies to the full time span. `WarnLevel` by default.
- `unbalanced_connector`: a `@connector` has a different number of flow variables and
  regular variables. `WarnLevel` by default.
- `affect_default_independent_variable`: an affect has no independent variable, so it
  uses `t_nounits`. `WarnLevel` by default.
- `substitute_skips_events`: `substitute(sys, rules)` changed a system that has events,
  but it does not change the events. `WarnLevel` by default. Set it with the `verbose`
  keyword of `substitute`.

# Verbosity stored on a system

Some diagnostics run where no `verbose` keyword is available: the checks in the
[`System`](@ref) constructor, `@connector`, callback affects, and `remake`. These use the
verbosity that is stored in the metadata of the system. The `verbose` keyword of `System`
sets it, and [`mtkcompile`](@ref) stores its own `verbose` value on the system that it
returns. Functions that have their own `verbose` keyword (problem constructors,
`mtkcompile`, `substitute`, ...) do not read the stored value.

# Sub-specifiers

- `initialization_verbosity`: the verbosity forwarded to the solve of the
  initialization `NonlinearProblem`/`SCCNonlinearProblem`. Accepts a `SciMLLogging`
  preset or a solver verbosity specifier (e.g. `NonlinearVerbosity`). Defaults to
  `Minimal()` (`None()` under the `None` preset), so inner linear-solver chatter stays
  off during initialization.

# Deprecated keyword equivalents

These boolean keywords remain accepted and, when explicitly passed, override the
corresponding toggles. They are deprecated in favor of `verbose`:

| Deprecated keyword | Toggles overridden |
|:--- |:--- |
| `warn_initialize_determined` | `singular_initialization`, `overdetermined_initialization`, `underdetermined_initialization` |
| `warn_cyclic_dependency` | `cyclic_dependency` |
| `warn_empty_op` (`linearization_function`) | `empty_operating_point` |

# Examples

```julia
using ModelingToolkit
using SciMLLogging: SciMLLogging, Silent, InfoLevel

mtkcompile(sys)                                                # Standard verbosity
mtkcompile(sys; verbose = false)                               # silence everything
ODEProblem(sys, op, tspan; verbose = SciMLLogging.None())      # quiet, MTK and solver
ODEProblem(sys, op, tspan; verbose = MTKVerbosity(initialization = Silent))
mtkcompile(sys; verbose = MTKVerbosity(compilation = InfoLevel))
```
""" MTKVerbosity

const DEFAULT_MTK_VERBOSE = MTKVerbosity()
const SILENT_MTK_VERBOSE = MTKVerbosity(None())

"""
    $(TYPEDEF)

Metadata key for the [`MTKVerbosity`](@ref) that is stored on a system. The `verbose`
keyword of `System` and [`mtkcompile`](@ref) set it. Diagnostics that have no `verbose`
keyword of their own read it with `_system_verbosity`.
"""
struct MTKVerbosityCtx end

# The verbosity stored in a metadata collection, or the default. `metadata` is the
# `ImmutableDict` of a system or a user-supplied collection of pairs. An `ImmutableDict`
# iterates the newest entry first, so the first match wins.
function _metadata_verbosity(metadata)
    for (k, v) in metadata
        k === MTKVerbosityCtx && return v::MTKVerbosity
    end
    return DEFAULT_MTK_VERBOSE
end

_system_verbosity(sys) = getmetadata(sys, MTKVerbosityCtx, DEFAULT_MTK_VERBOSE)::MTKVerbosity

@inline _process_verbose_param(verbose::MTKVerbosity) = verbose
@inline _process_verbose_param(preset::AbstractVerbosityPreset) = MTKVerbosity(preset)
@inline function _process_verbose_param(verbose::Bool)
    return verbose ? DEFAULT_MTK_VERBOSE : SILENT_MTK_VERBOSE
end

# Type-routed problem-constructor `verbose` (shared kwarg with the solver):
#   MTKVerbosity -> consumed by MTK; `filter_kwargs` drops it from `prob.kwargs`
#   preset       -> consumed (`MTKVerbosity(preset)`) AND forwarded to `prob.kwargs`
#   Bool         -> consumed (Standard/None) AND forwarded (long-standing behavior)
#   anything else (`nothing`, `DEVerbosity`, ...) -> MTK default; forwarded untouched
@inline function _route_problem_verbose(
        v::Union{MTKVerbosity, AbstractVerbosityPreset, Bool}
    )
    return _process_verbose_param(v)
end
@inline _route_problem_verbose(::Any) = DEFAULT_MTK_VERBOSE

# Whether a toggle can emit at all — used to guard work that only serves the message
# (e.g. substitution-cycle searches or printed reports).
@inline _toggle_enabled(::MTKVerbosity{false}, ::Symbol) = false
@inline function _toggle_enabled(verb::MTKVerbosity, toggle::Symbol)
    return verbosity_to_bool(getproperty(verb, toggle)::MessageLevel)
end

"""
    $(TYPEDSIGNATURES)

Apply a deprecated boolean keyword to `verb`: `nothing` (not explicitly passed) is a
no-op; `true` rebuilds the specifier with the toggles named in `pairs` forced to the
paired level (and enabled, so the override wins even over `MTKVerbosity(None())`);
`false` forces them to `Silent`.
"""
_override_toggle(verb::MTKVerbosity, ::Nothing, pairs::Pair...) = verb
function _override_toggle(
        verb::MTKVerbosity, flag::Bool, pairs::Pair{Symbol, MessageLevel}...
    )
    vals = map(fieldnames(typeof(verb))) do f
        i = findfirst(p -> first(p) === f, pairs)
        i === nothing ? getfield(verb, f) : (flag ? last(pairs[i]) : Silent)
    end
    E = flag ? true : verb isa MTKVerbosity{true}
    return MTKVerbosity{E}(vals...)
end

# The inner initialization problem inherits the outer cyclic-dependency setting: pass
# an explicit `true` override when the toggle is enabled, or `nothing` (no override).
function _cyclic_flag_for_initialization(verb::MTKVerbosity)
    return _toggle_enabled(verb, :cyclic_dependency) ? true : nothing
end
