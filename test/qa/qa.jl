using ModelingToolkit
using ModelingToolkitBase
using ModelingToolkitTearing
using StateSelection
using Aqua
using JET
using SciMLTesting
using Test

# ExplicitImports only sees an extension module once its trigger package is loaded, so
# every weakdep is loaded here to bring the extensions into the checked module set.
using FMIImport
using OrdinaryDiffEqBDF
using OrdinaryDiffEqDefault
using OrdinaryDiffEqRosenbrock

const MTK_EXTENSIONS = (
    :MTKFMIExt,
    :MTKOrdinaryDiffEqBDFExt,
    :MTKOrdinaryDiffEqDefaultExt,
    :MTKOrdinaryDiffEqRosenbrockExt,
)

# ExplicitImports silently skips an extension that fails to load, so assert the
# extension modules actually exist rather than trusting a green run_qa.
@testset "Extensions loaded" begin
    for ext in MTK_EXTENSIONS
        @test Base.get_extension(ModelingToolkit, ext) !== nothing
    end
end

# The exact set of externally-owned bindings ModelingToolkit deliberately exposes as part
# of its own public API. ModelingToolkit is a facade: `@reexport using ModelingToolkitBase`
# pulls up the modelling layer, which in turn `@reexport`s Symbolics (and, through it,
# SymbolicUtils and TermInterface), and the SciMLBase problem types are what every
# documented workflow constructs. Without this list `public_reexports` reports all 487 of
# them, because a blanket reexport cannot distinguish "we mean to expose this" from "a
# dependency happened to add a public binding".
#
# Pinning the set here is the point: it is version-controlled, so a dependency adding a new
# public name fails this test instead of silently widening ModelingToolkit's API, and every
# addition becomes a deliberate decision. Removing an entry is a breaking change.
const REEXPORTED_API = (
    # ModelingToolkitBase: the lower half of this library. ModelingToolkit is a facade
    # over it (`@reexport using ModelingToolkitBase` plus `@import_mtkbase`), so its entire
    # public API is deliberately part of ModelingToolkit's.
    Symbol("@brownian"), Symbol("@brownians"), Symbol("@component"), Symbol("@connector"),
    Symbol("@constants"), Symbol("@discretes"), Symbol("@independent_variables"),
    Symbol("@mtkbuild"), Symbol("@mtkcompile"), Symbol("@mtkcomplete"), Symbol("@named"),
    Symbol("@namespace"), Symbol("@nonamespace"), Symbol("@parameters"),
    Symbol("@poissonians"), :AbstractCollocation, :AbstractSystem, :add_accumulations,
    :alg_equations, :AnalysisPoint, :analytically_integrated, :apply_to_variables, :asdigraph,
    :asgraph, :assertions, :AssignmentAffect, :bindings, :Both, :bound_inputs, :bound_outputs,
    :bound_parameters, :brownians, :calculate_control_jacobian, :calculate_cost_gradient,
    :calculate_cost_hessian, :calculate_hessian, :calculate_jacobian, :calculate_massmatrix,
    :calculate_tgrad, :CasADiCollocation, :CasADiDynamicOptProblem,
    :change_independent_variable, :change_of_variables, :check_mutable_cache,
    :check_symbolic_ad_allowed, :CheckAll, :CheckComponents, :CheckNone, :CheckUnits,
    :collect_scoped_vars!, :collect_var_to_name!, :collect_vars!, :CompilerOptions, :complete,
    :compose, :connect, :Connection, :constraints, :continuous_events,
    :continuous_events_toplevel, :convert_bindings_for_time_independent_system,
    :convert_system_indepvar, :cost, :D, :D_nounits, :debug_system, :diff_equations,
    :discrete_events, :discrete_events_toplevel, :DiscreteSystem, :domain_connect,
    :DynamicOptSolution, :eqeq_dependencies, :eqtype_supports_collect_vars, :Equality,
    :equation_dependencies, :equations, :equations_toplevel, :EvalAt, :expand_connections,
    :extend, :flatten, :Flow, :fractional_to_ordinary, :full_equations,
    :generate_control_jacobian, :generate_cost, :generate_cost_gradient,
    :generate_cost_hessian, :generate_custom_function, :generate_diffusion_function,
    :generate_initializesystem, :generate_jacobian, :generate_rhs, :generate_tgrad,
    :generate_W, :get_alg_eqs, :get_analytically_integrated, :get_assertions, :get_bcs,
    :get_bindings, :get_brownians, :get_connector_type, :get_consolidate, :get_constraints,
    :get_continuous_events, :get_costs, :get_description, :get_diff_eqs, :get_discrete_events,
    :get_domain, :get_dvs, :get_eqs, :get_guesses, :get_gui_metadata,
    :get_ignored_connections, :get_index_cache, :get_initial_conditions,
    :get_initialization_eqs, :get_initializesystem, :get_inputs, :get_irreducibles,
    :get_irstructure_tlv, :get_is_dde, :get_is_discrete, :get_is_initializesystem,
    :get_isscheduled, :get_iv, :get_ivs, :get_jumps, :get_maybe_zeros, :get_metadata,
    :get_name, :get_noise_eqs, :get_observed, :get_outputs, :get_parameter_bindings_graph,
    :get_parent, :get_poissonians, :get_preface, :get_ps, :get_schedule,
    :get_state_priorities, :get_systems, :get_tag, :get_tearing_state, :get_tspan,
    :get_tstops, :get_unknowns, :get_var_to_name, :get_w, :getbounds, :getconnect,
    :getdefault, :getdescription, :getdist, :getguess, :getmisc, :getnominal, :getunit,
    :Girsanov_transform, :GlobalScope, :guesses, :has_alg_eqs, :has_alg_equations,
    :has_analytically_integrated, :has_assertions, :has_bcs, :has_bindings, :has_brownians,
    :has_connector_type, :has_consolidate, :has_constraints, :has_continuous_events,
    :has_costs, :has_description, :has_diff_eqs, :has_diff_equations, :has_discrete_events,
    :has_domain, :has_dvs, :has_eqs, :has_guesses, :has_gui_metadata,
    :has_ignored_connections, :has_index_cache, :has_initial_conditions,
    :has_initialization_eqs, :has_initializesystem, :has_inputs, :has_irreducibles,
    :has_irstructure_tlv, :has_is_dde, :has_is_discrete, :has_is_initializesystem,
    :has_isscheduled, :has_iv, :has_ivs, :has_jumps, :has_maybe_zeros, :has_metadata,
    :has_name, :has_noise_eqs, :has_observed, :has_outputs, :has_parameter_bindings_graph,
    :has_parent, :has_poissonians, :has_preface, :has_ps, :has_schedule,
    :has_state_priorities, :has_systems, :has_tag, :has_tearing_state, :has_tspan,
    :has_tstops, :has_unknowns, :has_var_to_name, :hasbounds, :hasconnect, :hasdefault,
    :hasdescription, :hasdist, :hasguess, :hasmisc, :hasnominal, :hasunit, :hierarchy, :Hold,
    :homotopy, :HomotopyContinuationProblem, :ImperativeAffect, :ImplicitDiscreteSystem,
    :independent_variable, :independent_variables, :InfiniteOptCollocation,
    :InfiniteOptDynamicOptProblem, :Initial, :initial_conditions, :initialization_equations,
    :InitializationProblem, :inputs, :instream, :irreducibles, :is_alg_equation, :is_bound,
    :is_diff_equation, :iscomplete, :isdisturbance, :isinitial, :isinput, :isirreducible,
    :isoutput, :isparameter, :istunable, :JuMPCollocation, :JuMPDynamicOptProblem, :jumps,
    :JumpSystem, :linear_fractional_to_ordinary, :liouville_transform, :LocalScope,
    :maybe_zeros, :MiscSystemData, :MissingGuessValue, :ModelingToolkitBase,
    :modelingtoolkitize, :modified_unknowns!, :mtkcompile, :MTKParameters,
    :MTKVariableTypeCtx, :namespace_equations, :noise_to_brownians, :NonlinearSystem,
    :observables, :observed, :ODESystem, :open_loop, :OptimizationSystem, :outputs,
    :parameters, :parameters_toplevel, :ParentScope, :PDESystem, :Pre, :ProblemTypeCtx,
    :PyomoCollocation, :PyomoDynamicOptProblem, :renamespace, :reorder_dimension_by_tunables,
    :reorder_dimension_by_tunables!, :respecialize, :Sample, :SampleTime, :SDESystem,
    :set_defaults, :setdefault, :setguess, :setnominal, :Shift, :ShiftIndex,
    :should_invalidate_mutable_cache_entry, :state_priorities, :state_priority,
    :stochastic_integral_transform, :store_to_mutable_cache!, :Stream, :structural_simplify,
    :subset_tunables, :SymbolicADDisallowed, :SymbolicContinuousCallback,
    :SymbolicDiscreteCallback, :SymbolicMassActionJump, :SymScope, :System, :t, :t_nounits,
    :tobrownian, :toggle_namespacing, :toparam, :tunable_parameters, :unbound_inputs,
    :unbound_outputs, :unknowns, :unknowns_toplevel, :variable_dependencies, :VariableBounds,
    :VariableConnectType, :VariableDescription, :VariableInput, :VariableIrreducible,
    :VariableMisc, :VariableOutput, :VariableStatePriority, :VariableType, :VariableUnit,
    :varmap_to_vars, :varvar_dependencies,
    # Symbolics: the symbolic algebra layer users write models in. Reached through
    # `@reexport using Symbolics` in ModelingToolkitBase.
    Symbol("@derivative_rule"), Symbol("@derivatives"), Symbol("@register_array_symbolic"),
    Symbol("@register_derivative"), Symbol("@register_discontinuity"),
    Symbol("@register_inverse"), Symbol("@register_symbolic"), Symbol("@symbolic_wrap"),
    Symbol("@symstruct"), Symbol("@variables"), Symbol("@wrapped"), :approximation_function,
    :build_function, :Differential, :Equation, :expand_derivatives, :factors, :gather_factor,
    :get_variables, :groebner_basis, :has_inverse, :has_left_inverse, :has_right_inverse,
    :Inequality, :infimum, :Integral, :inverse, :inverse_laplace, :is_derivative,
    :is_groebner_basis, :laplace, :laplace_solve_ode, :left_continuous_function,
    :left_inverse, :limit, :majorization_function, :minorization_function, :Num,
    :parse_expr_to_symbolic, :partial_frac_decomposition, :polynomial_coeffs,
    :right_continuous_function, :right_inverse, :rootfunction, :semilinear_form,
    :semipolynomial_form, :semiquadratic_form, :series, :solve_for, :solve_linear_ode_system,
    :solve_symbolic_IVP, :substitute_in_deriv, :substitute_in_deriv_and_depvar, :supremum,
    :symbolic_linear_solve, :symbolic_solve, :symbolic_solve_ode, :SymbolicLinearODE,
    :Symbolics, :symbolics_to_sympy, :symbolics_to_sympy_pythoncall,
    :SymbolicsSparsityDetector, :sympy_algebraic_solve, :sympy_integrate, :sympy_limit,
    :sympy_linear_solve, :sympy_ode_solve, :sympy_pythoncall_algebraic_solve,
    :sympy_pythoncall_integrate, :sympy_pythoncall_limit, :sympy_pythoncall_linear_solve,
    :sympy_pythoncall_ode_solve, :sympy_pythoncall_simplify, :sympy_pythoncall_to_symbolics,
    :sympy_simplify, :sympy_to_symbolics, :SymStruct, :taylor, :taylor_coeff, :terms,
    :tosymbol, :≲, :≳,
    # SymbolicUtils (including its `Rewriters` and `Code` submodules): re-exported
    # transitively by Symbolics.
    Symbol("@acrule"), Symbol("@arrayop"), Symbol("@makearray"), Symbol("@rule"),
    Symbol("@syms"), :BS, :expand, :flatten_fractions, :get_canonical_expr, :get_reachability,
    :getmetadata, :hasmetadata, :ifelse_branching, :ifelse_eager, :IRStructure, :istree,
    :populate_ir!, :print_ir, :quick_cancel, :Rewriters, :RuleSet, :SafeReal, :setmetadata,
    :simplify, :simplify_fractions, :substitute, :SymbolicUtils, :SymReal, :Term, :term,
    :toexpr, :TreeReal, :unwrap_const, :vartype,
    # TermInterface: the term-manipulation interface, re-exported transitively by
    # SymbolicUtils.
    :arguments, :iscall, :operation, :sorted_arguments,
    # SciMLBase: the problem/function types ModelingToolkit constructs. Users call
    # `ODEProblem(sys, ...)` unqualified after `using ModelingToolkit`.
    :AbstractDynamicOptProblem, :AbstractNonlinearProblem, :Clock, :DiscreteFunction,
    :DiscreteProblem, :ImplicitDiscreteFunction, :ImplicitDiscreteProblem,
    :IntervalNonlinearFunction, :IntervalNonlinearProblem, :NonlinearFunction,
    :NonlinearProblem, :ODEFunction, :ODEProblem, :OptimizationProblem, :SDEFunction,
    :SDEProblem, :SolverStepClock, :SteadyStateProblem, :TimeDomain,
    # CommonSolve: `solve` is the entry point of every documented workflow.
    :solve,
    # JumpProcesses: `JumpProblem` is constructed from a `System` with jumps.
    :JumpProblem,
    # BipartiteGraphs: `BipartiteGraph` is re-exported by ModelingToolkitBase for the
    # dependency-graph API (`equation_dependencies` and friends return one).
    :BipartiteGraph,
    # ModelingToolkitTearing: `TearingState` is explicitly exported by
    # ModelingToolkit (see `src/ModelingToolkit.jl`).
    :TearingState,
    # StateSelection: `find_solvables!` is an accidental leak -- it is a StateSelection
    # internal that `StructuralTransformations` `export`s but only ever calls qualified. It is
    # listed here so the audit is exhaustive; dropping it is a breaking change and belongs in
    # the next major release, not in a QA fix.
    :find_solvables!,
    # UnPack: `@unpack`/`@pack!`, re-exported by ModelingToolkitBase.
    Symbol("@pack!"), Symbol("@unpack"), :UnPack,
)

# Public names that reach ModelingToolkit's API surface only through
# `@reexport using Symbolics` in ModelingToolkitBase. They are owned (and undocumented) by
# Symbolics/SymbolicUtils, so ModelingToolkit is not the right place to document them.
const SYMBOLICS_OWNED_REEXPORTS = (
    Symbol("@symbolic_wrap"),
    Symbol("@wrapped"),
    :RuleSet,
    :get_canonical_expr,
    :infimum,
    :is_derivative,
    :istree,
    :solve_for,
    :supremum,
)

# Names the extensions must reach for which no public spelling exists.
const NONPUBLIC_QUALIFIED_ACCESSES = (
    # ModelingToolkitBase: the FMI extension's variable metadata helpers.
    :default_toterm,
    :setinput,
    :setoutput,
    # ModelingToolkit's own stub that MTKFMIExt implements; it is documented API but has
    # not been declared `public`, and an extension has no way to add a method unqualified.
    :FMIComponent,
    # Base: `Base.RefValue` has no public spelling.
    :RefValue,
)

# ModelingToolkit is the upper half of ModelingToolkitBase: the two are one library split
# across a monorepo boundary, with the compilation and problem-construction layer living
# here and dispatching on the system and operator types defined below it. Methods added to
# those types are extensions of the pair's own API, so Aqua should not count them as piracy.
const MTKBASE_OWNED_TYPES = (
    ModelingToolkitBase.Hold,
    ModelingToolkitBase.Sample,
    ModelingToolkitBase.SampleTime,
    ModelingToolkitBase.ShiftIndex,
    ModelingToolkitBase.System,
)

# ModelingToolkit is also where the structural-transformation stack is wired to
# ModelingToolkitBase's hooks, which means forwarding a handful of MTKBase functions to the
# StateSelection/ModelingToolkitTearing implementations on their own types. `TearingState`
# is not yet public in ModelingToolkitTearing, so it has to be reached by qualified access.
const STRUCTURAL_TYPES = (
    ModelingToolkitTearing.TearingState,
    StateSelection.DiffGraph,
)

run_qa(
    ModelingToolkit;
    Aqua = Aqua,
    JET = JET,
    jet = true,
    jet_kwargs = (; target_defined_modules = true),
    aqua_kwargs = (;
        piracies = (; treat_as_own = (MTKBASE_OWNED_TYPES..., STRUCTURAL_TYPES...)),
    ),
    ei_kwargs = (;
        all_qualified_accesses_are_public = (; ignore = NONPUBLIC_QUALIFIED_ACCESSES),
    ),
    api_docs_kwargs = (; ignore = SYMBOLICS_OWNED_REEXPORTS),
    reexports_allow = REEXPORTED_API,
)
