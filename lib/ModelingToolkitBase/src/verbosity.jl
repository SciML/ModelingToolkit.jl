@verbosity_specifier MTKVerbosity begin
    sub_specifiers = (
        :initialization_verbosity, :de_verbosity, :nonlinear_verbosity,
        :linear_verbosity, :optimization_verbosity, :bvp_verbosity,
    )

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
            de_verbosity = None(),
            nonlinear_verbosity = None(),
            linear_verbosity = None(),
            optimization_verbosity = None(),
            bvp_verbosity = None(),
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
            de_verbosity = Minimal(),
            nonlinear_verbosity = Minimal(),
            linear_verbosity = Minimal(),
            optimization_verbosity = Minimal(),
            bvp_verbosity = Minimal(),
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
            de_verbosity = Standard(),
            nonlinear_verbosity = Standard(),
            linear_verbosity = Standard(),
            optimization_verbosity = Standard(),
            bvp_verbosity = Standard(),
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
            de_verbosity = Detailed(),
            nonlinear_verbosity = Detailed(),
            linear_verbosity = Detailed(),
            optimization_verbosity = Detailed(),
            bvp_verbosity = Detailed(),
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
            de_verbosity = All(),
            nonlinear_verbosity = All(),
            linear_verbosity = All(),
            optimization_verbosity = All(),
            bvp_verbosity = All(),
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

- an `MTKVerbosity` sets ModelingToolkit's verbosity, and its solver sub-specifier for
  the problem type (see below) is placed in `prob.kwargs`. For example, an `ODEProblem`
  gets `verbose = v.de_verbosity`;
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
  regular variables. `mtkcompile` reports it once for each connector type. The
  `verbose` of `mtkcompile` and the verbosity stored on the connector must both enable
  it. `WarnLevel` by default.
- `affect_default_independent_variable`: an affect has no independent variable, so it
  uses `t_nounits`. `WarnLevel` by default.
- `substitute_skips_events`: `substitute(sys, rules)` changed a system that has events,
  but it does not change the events. `WarnLevel` by default. Set it with the `verbose`
  keyword of `substitute`.

# Verbosity stored on a system

Some diagnostics run where no `verbose` keyword is available: the checks in the
[`System`](@ref) constructor, callback affects, and `remake`. These use the
verbosity that is stored in the metadata of the system. The `verbose` keyword of `System`
sets it, and [`mtkcompile`](@ref) stores its own `verbose` value on the system that it
returns. If `mtkcompile` gets no `verbose`, it uses the value stored on its input system.
Other functions that have their own `verbose` keyword (problem constructors,
`substitute`, ...) do not read the stored value.

# Sub-specifiers

- `initialization_verbosity`: the verbosity forwarded to the solve of the
  initialization `NonlinearProblem`/`SCCNonlinearProblem`. Accepts a `SciMLLogging`
  preset or a solver verbosity specifier (e.g. `NonlinearVerbosity`). Defaults to
  `Minimal()` (`None()` under the `None` preset), so inner linear-solver chatter stays
  off during initialization. A preset also applies to a linear initialization problem;
  a `NonlinearVerbosity` does not, and that problem uses `linear_verbosity` instead.
  The ModelingToolkit diagnostics of the initialization problem use the toggles of the
  outer `MTKVerbosity`.

The other sub-specifiers hold the verbosity of the solver for each problem family. A
problem constructor that gets an `MTKVerbosity` forwards the matching field to `solve`.
Each field accepts a `SciMLLogging` preset or the solver's own specifier, and each preset
of `MTKVerbosity` sets them to the same preset (for example, `MTKVerbosity(None())` has
`de_verbosity = None()`).

| Field | Problem types | Solver specifier |
|:--- |:--- |:--- |
| `de_verbosity` | `ODEProblem`, `DAEProblem`, `SDEProblem`, `DDEProblem`, `SDDEProblem`, `DiscreteProblem`, `ImplicitDiscreteProblem`, `JumpProblem` | `DEVerbosity` (DiffEqBase) |
| `nonlinear_verbosity` | `NonlinearProblem`, `NonlinearLeastSquaresProblem`, `IntervalNonlinearProblem`, `SteadyStateProblem`, `HomotopyProblem` | `NonlinearVerbosity` (NonlinearSolveBase) |
| `linear_verbosity` | `LinearProblem` | `LinearVerbosity` (LinearSolve) |
| `optimization_verbosity` | `OptimizationProblem` | `OptimizationVerbosity` (OptimizationBase) |
| `bvp_verbosity` | `BVProblem` | `BVPVerbosity` (BoundaryValueDiffEqCore) |

A `SteadyStateProblem` solved with `DynamicSS` does not give `nonlinear_verbosity` to the
ODE solve; use the `odesolve_kwargs` of `DynamicSS` for that solve. The dynamic
optimization problems (for example `JuMPDynamicOptProblem`) store no solver verbosity,
because their `solve` methods take their own `verbose::Bool`.

The solver specifiers have their own sub-specifiers. For example, the linear solves
inside an ODE solve are controlled by the `linear_verbosity` of the `DEVerbosity`, not by
the `linear_verbosity` field of `MTKVerbosity`:

```julia
MTKVerbosity(de_verbosity = DEVerbosity(linear_verbosity = SciMLLogging.None()))
```

# Deprecated keyword equivalents

These boolean keywords are deprecated in favor of `verbose`. When explicitly passed,
they emit a deprecation warning and override the corresponding toggles:

| Deprecated keyword | Toggles overridden |
|:--- |:--- |
| `warn_initialize_determined` | `singular_initialization`, `overdetermined_initialization`, `underdetermined_initialization` |
| `warn_cyclic_dependency` | `cyclic_dependency` |
| `warn_empty_op` (`linearization_function`) | `empty_operating_point` |
| `print_underconstrained_variables` (`mtkcompile`) | `underconstrained_variables` (`true` sets `InfoLevel`) |

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
keyword of their own read it with `getmetadata(sys, MTKVerbosityCtx, MTKVerbosity())`.
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

# `sys` with `verbose` stored. No new copy of `sys` is made if it already stores `verbose`.
function _set_system_verbosity(sys, verbose::MTKVerbosity)
    getmetadata(sys, MTKVerbosityCtx, nothing) === verbose && return sys
    return setmetadata(sys, MTKVerbosityCtx, verbose)
end

# `MTKVerbosity(verbose)` accepts every value that a `verbose = ...` keyword of
# ModelingToolkit documents: an `MTKVerbosity`, a preset, or a `Bool`.
MTKVerbosity(verbose::MTKVerbosity) = verbose
MTKVerbosity(verbose::Bool) = verbose ? DEFAULT_MTK_VERBOSE : SILENT_MTK_VERBOSE

@inline _process_verbose_param(verbose::MTKVerbosity) = verbose
@inline _process_verbose_param(preset::AbstractVerbosityPreset) = MTKVerbosity(preset)
@inline function _process_verbose_param(verbose::Bool)
    return verbose ? DEFAULT_MTK_VERBOSE : SILENT_MTK_VERBOSE
end

# Type-routed problem-constructor `verbose` (shared kwarg with the solver):
#   MTKVerbosity -> used by MTK; `_solver_kwargs` puts its solver sub-specifier for the
#                   problem type (e.g. `de_verbosity`) in `prob.kwargs`
#   preset       -> used as `MTKVerbosity(preset)` AND forwarded to `prob.kwargs`
#   Bool         -> used as Standard/None AND forwarded (long-standing behavior)
#   anything else (`nothing`, `DEVerbosity`, ...) -> MTK default; forwarded untouched
@inline function _route_problem_verbose(
        v::Union{MTKVerbosity, AbstractVerbosityPreset, Bool}
    )
    return _process_verbose_param(v)
end
@inline _route_problem_verbose(::Any) = DEFAULT_MTK_VERBOSE

# The toggles that the deprecated `warn_initialize_determined` keyword overrides.
const _INITIALIZE_DETERMINED_TOGGLES = (
    :singular_initialization => WarnLevel,
    :overdetermined_initialization => WarnLevel,
    :underdetermined_initialization => WarnLevel,
)

# The `MTKVerbosity` of a problem constructor with keywords `kwargs`: `verbose` with the
# deprecated boolean keywords applied, as in `SciMLProblemOptions`. Use it where the
# `SciMLProblemOptions` of the problem is not available. `SciMLProblemOptions` emits the
# deprecation warnings, so this function does not.
function _problem_verbosity(kwargs)
    verbosity = _route_problem_verbose(get(kwargs, :verbose, nothing))
    verbosity = _override_toggle(
        verbosity, get(kwargs, :warn_initialize_determined, nothing),
        _INITIALIZE_DETERMINED_TOGGLES...
    )
    return _override_toggle(
        verbosity, get(kwargs, :warn_cyclic_dependency, nothing),
        :cyclic_dependency => WarnLevel
    )
end

# The `MTKVerbosity` sub-specifier that holds the solver verbosity for problem type `T`.
# `AbstractBVProblem <: AbstractODEProblem`, so it needs its own method.
_solver_verbosity_field(::Type) = nothing
_solver_verbosity_field(::Type{<:SciMLBase.AbstractDEProblem}) = :de_verbosity
_solver_verbosity_field(::Type{<:SciMLBase.AbstractBVProblem}) = :bvp_verbosity
# The dynamic optimization `solve` methods take their own `verbose::Bool` and do not read
# `prob.kwargs`, so an `MTKVerbosity` is not stored.
_solver_verbosity_field(::Type{<:SciMLBase.AbstractDynamicOptProblem}) = nothing
_solver_verbosity_field(::Type{<:SciMLBase.AbstractNonlinearProblem}) = :nonlinear_verbosity
_solver_verbosity_field(::Type{<:SciMLBase.AbstractLinearProblem}) = :linear_verbosity
function _solver_verbosity_field(::Type{<:SciMLBase.AbstractOptimizationProblem})
    return :optimization_verbosity
end

# Replace an `MTKVerbosity` in the keywords of a problem of type `T` by the solver
# sub-specifier for `T`, so that `solve` gets a value it accepts. Other values of
# `verbose` are not changed.
function _solver_kwargs(::Type{T}, kwargs) where {T}
    verbose = get(kwargs, :verbose, nothing)
    verbose isa MTKVerbosity || return kwargs
    nt = NamedTuple(kwargs)
    rest = NamedTuple{filter(!=(:verbose), keys(nt))}(nt)
    field = _solver_verbosity_field(T)
    field === nothing && return rest
    return merge(rest, (; verbose = getfield(verbose, field)))
end

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

# `_override_toggle` for the deprecated boolean keyword `kw` of `caller`, with a
# deprecation warning when the keyword is explicitly passed.
function _deprecated_override_toggle(
        verb::MTKVerbosity, ::Nothing, kw::Symbol, caller::Symbol, pairs::Pair...
    )
    return verb
end
function _deprecated_override_toggle(
        verb::MTKVerbosity, flag::Bool, kw::Symbol, caller::Symbol,
        pairs::Pair{Symbol, MessageLevel}...
    )
    toggles = join(("$(first(p)) = ..." for p in pairs), ", ")
    Base.depwarn(
        "The `$kw` keyword is deprecated, use `verbose = MTKVerbosity($toggles)` instead.",
        caller
    )
    return _override_toggle(verb, flag, pairs...)
end

# `verb` with `solver` as its `initialization_verbosity`. `nothing` and an `MTKVerbosity`
# do not change `verb`.
function _with_initialization_solver_verbosity(
        verb::MTKVerbosity, ::Union{Nothing, MTKVerbosity}
    )
    return verb
end
function _with_initialization_solver_verbosity(verb::MTKVerbosity{E}, solver) where {E}
    vals = map(fieldnames(typeof(verb))) do f
        f === :initialization_verbosity ? solver : getfield(verb, f)
    end
    return MTKVerbosity{E}(vals...)
end

# The `verbose` of the inner initialization problem: the toggles of `verb`, with its
# `initialization_verbosity` as the solver verbosity of the nonlinear problem types. A
# preset or `Bool` also becomes the `linear_verbosity`; a `NonlinearVerbosity` does not,
# because `LinearProblem` does not accept it.
function _initialization_problem_verbosity(verb::MTKVerbosity{E}) where {E}
    init = verb.initialization_verbosity
    linear = init isa Union{AbstractVerbosityPreset, Bool} ? init : verb.linear_verbosity
    vals = map(fieldnames(typeof(verb))) do f
        f === :nonlinear_verbosity ? init :
            f === :linear_verbosity ? linear : getfield(verb, f)
    end
    return MTKVerbosity{E}(vals...)
end
