using ModelingToolkitBase
using Aqua
using JET
using SciMLTesting

# This is loaded only by MTKBifurcationKitExt, while remaining a hard dependency so
# loading BifurcationKit alone still activates the extension.
const STALE_EXTENSION_DEPENDENCIES = [:SimpleNonlinearSolve]

# ModelingToolkitBase supplies Symbolics expression-protocol methods for equations and
# JumpProcesses jump types. These methods intentionally extend the corresponding
# SymbolicUtils generic functions, while all other piracies remain checked.
const INTENTIONAL_EXTERNAL_GENERIC_EXTENSIONS = (
    ModelingToolkitBase.SymbolicUtils.search_variables!,
    ModelingToolkitBase.toexpr,
)

# Moshi's `@data` expands to a submodule whose source ExplicitImports cannot walk.
# SciMLTesting already allows the equivalent `EnumX.@enumx` modules automatically.
const GENERATED_ADT_MODULES = (
    ModelingToolkitBase.MissingGuessValue,
    ModelingToolkitBase.StructuralHint,
)

# Externally-owned names ModelingToolkitBase imports for which the owning package offers
# no public spelling.
#
#   Symbolics/SymbolicUtils. ModelingToolkitBase is built directly on the Symbolics term
#   representation: the `BasicSymbolic` implementation types (`SArgsT`, `SConst`, `SSym`,
#   `STerm`, `SymbolicT`, `VartypeT`, `BSImpl`, `Term`, `Operator`, `CallAndWrap`), the
#   sparsity/linearity analyses and the naming helpers all live below the documented
#   Symbolics surface.
#
#   StaticArraysCore, FunctionWrappers, OffsetArrays. `StaticArray`, `StaticVector` and
#   `similar_type` are declared by StaticArraysCore but exported only by StaticArrays, so
#   the owner and the public spelling disagree; `FunctionWrapper` and `Origin` have no
#   public spelling at all.
const NONPUBLIC_EXPLICIT_IMPORTS = (
    # Symbolics
    :CallAndWrap, Symbol("@derivatives"), :fixpoint_sub, :hasderiv, :isaffine, :islinear,
    :jacobian_sparsity, :lower_varname, :recursive_hasoperator, :rename, :SArgsT, :SConst,
    :setdefaultval, :SSym, :STerm, :SymbolicT, :var_from_nested_derivative, :VartypeT,
    # SymbolicUtils
    :_iszero, :BSImpl, :Operator, :Term,
    # StaticArraysCore
    :similar_type, :StaticArray, :StaticVector,
    # FunctionWrappers
    :FunctionWrapper,
    # OffsetArrays
    :Origin,
)

# Externally-owned names reached by qualified access for which the owning package offers no
# public spelling. Beyond the Symbolics/SymbolicUtils internals covered above, these are:
# the SciMLBase callback/keyword defaults every problem constructor has to reproduce, the
# `Moshi.Data`/`Moshi.Match` macro entry points, the Base internals with no public spelling
# (`@__doc__`, `@ntuple`, `literal_pow`, `ReshapedArray`, ...), and the AD hooks
# (`ForwardDiff.Dual`, `EnzymeCore.EnzymeRules`) the generated code dispatches on.
const NONPUBLIC_QUALIFIED_ACCESSES = (
    # Symbolics
    :CallAndWrap, :canonical_form, :COMMON_ZERO, :executediff, :fixpoint_sub,
    :FixpointSubstituter, :FPSubFilterer, :getdefaultval, :getmetadata_maybe_indexed,
    :hessian_sparsity, :hide_lhs, :ia_solve, :isarraysymbolic, :iswrapped,
    :jacobian_sparsity, :leq, :parse_vars, :rename, :RootsOf, :SArgsT, :SConst, :SSym,
    :STerm, :warn_load_latexify,
    # SymbolicUtils
    :_iszero, :ACDict, :AddMulVariant, :array_literal, :ArrayMaker, :as_linear_idx,
    Symbol("@cache"), :Code, :Const, :default_is_atomic, :default_substitute_filter, :Fill,
    :hashcons, :idxs_for_arrayop, :IRStructureSearchBuffer, :IRSubstituter, :is_array_shape,
    :is_called_function_symbolic, :is_function_symbolic, :isbinop, :isequal_somescalar,
    :Operator, :promote_symtype, :query, :RecursiveDFS, :RegionsT, :search_variables,
    :search_variables!, :SmallV, :stable_eachindex, :StableIndex, :StableIndices,
    :subset_ir, :Substituter, :TypeT,
    # SymbolicUtils.Code
    :create_array,
    # SymbolicIndexingInterface
    :supports_tuple_observed,
    # SciMLBase
    :allowedkeywords, :anyeltypedual, :diagnose_symbolic_instability,
    :DISCRETE_INPLACE_DEFAULT, :FINALIZE_DEFAULT, :INITIALIZE_DEFAULT,
    # JumpProcesses
    :NullAggregator,
    # BipartiteGraphs
    :print_hyperedge_hint,
    # AbstractTrees
    :printnode,
    # BandedMatrices
    :_BandedMatrix,
    # EnzymeCore
    :EnzymeRules,
    # EnzymeCore.EnzymeRules
    :inactive_noinl, :inactive_type,
    # ForwardDiff
    :Dual, :valtype,
    # Moshi
    :Data, :Derive, :Match,
    # Setfield
    :PropertyLens,
    # StaticArraysCore
    :similar_type,
    # UnPack
    :unpack,
    # Base
    Symbol("@__doc__"), :Callable, :Cartesian, :deepcopy_internal, :Experimental,
    :HasEltype, :JLOptions, :literal_pow, Symbol("@nospecializeinfer"), Symbol("@ntuple"),
    :ReshapedArray,
    # Base.Cartesian
    :poplinenum,
    # Base.Experimental
    Symbol("@compiler_options"),
    # Base.Meta
    :ParseError,
)

# The exact set of externally-owned bindings ModelingToolkitBase deliberately exposes as
# part of its own public API. `@reexport using Symbolics` (and, through it, SymbolicUtils
# and TermInterface) pulls up the whole symbolic layer, `@reexport using UnPack` the
# destructuring macros, and the SciMLBase problem types are what every documented workflow
# constructs.
#
# Pinning the set here is the point: it is version-controlled, so a dependency adding a new
# public name fails this test instead of silently widening ModelingToolkitBase's API, and
# every addition becomes a deliberate decision. Removing an entry is a breaking change.
const REEXPORTED_API = (
    # Symbolics
    :approximation_function, :build_function, Symbol("@derivative_rule"),
    Symbol("@derivatives"), :Differential, :Equation, :expand_derivatives, :factors,
    :gather_factor, :get_variables, :groebner_basis, :has_inverse, :has_left_inverse,
    :has_right_inverse, :Inequality, :infimum, :Integral, :inverse, :inverse_laplace,
    :is_derivative, :is_groebner_basis, :laplace, :laplace_solve_ode,
    :left_continuous_function, :left_inverse, :limit, :majorization_function,
    :minorization_function, :Num, :parse_expr_to_symbolic, :partial_frac_decomposition,
    :polynomial_coeffs, Symbol("@register_array_symbolic"), Symbol("@register_derivative"),
    Symbol("@register_discontinuity"), Symbol("@register_inverse"),
    Symbol("@register_symbolic"), :right_continuous_function, :right_inverse, :rootfunction,
    :semilinear_form, :semipolynomial_form, :semiquadratic_form, :series, :solve_for,
    :solve_linear_ode_system, :solve_symbolic_IVP, :substitute_in_deriv,
    :substitute_in_deriv_and_depvar, :supremum, :symbolic_linear_solve, :symbolic_solve,
    :symbolic_solve_ode, Symbol("@symbolic_wrap"), :SymbolicLinearODE, :Symbolics,
    :symbolics_to_sympy, :symbolics_to_sympy_pythoncall, :SymbolicsSparsityDetector,
    :sympy_algebraic_solve, :sympy_integrate, :sympy_limit, :sympy_linear_solve,
    :sympy_ode_solve, :sympy_pythoncall_algebraic_solve, :sympy_pythoncall_integrate,
    :sympy_pythoncall_limit, :sympy_pythoncall_linear_solve, :sympy_pythoncall_ode_solve,
    :sympy_pythoncall_simplify, :sympy_pythoncall_to_symbolics, :sympy_simplify,
    :sympy_to_symbolics, Symbol("@symstruct"), :SymStruct, :taylor, :taylor_coeff, :terms,
    :tosymbol, Symbol("@variables"), Symbol("@wrapped"), :≲, :≳,
    # SymbolicUtils
    Symbol("@acrule"), Symbol("@arrayop"), :BS, :expand, :flatten_fractions,
    :get_reachability, :getmetadata, :hasmetadata, :ifelse_branching, :ifelse_eager,
    :IRStructure, :istree, Symbol("@makearray"), :populate_ir!, :print_ir, :quick_cancel,
    :Rewriters, Symbol("@rule"), :RuleSet, :SafeReal, :scalarize, :setmetadata, :shape,
    :simplify, :simplify_fractions, :substitute, :SymbolicUtils, :SymReal, Symbol("@syms"),
    :Term, :term, :TreeReal, :Unknown, :unwrap, :unwrap_const, :vartype,
    # SymbolicUtils.Code
    :toexpr,
    # TermInterface
    :arguments, :iscall, :operation, :sorted_arguments,
    # SciMLBase
    :AbstractDynamicOptProblem, :AbstractNonlinearProblem, :Clock, :DiscreteFunction,
    :DiscreteProblem, :ImplicitDiscreteFunction, :ImplicitDiscreteProblem,
    :IntervalNonlinearFunction, :IntervalNonlinearProblem, :NonlinearFunction,
    :NonlinearProblem, :ODEFunction, :ODEProblem, :OptimizationProblem, :SDEFunction,
    :SDEProblem, :SolverStepClock, :SteadyStateProblem, :TimeDomain,
    # CommonSolve
    :solve,
    # JumpProcesses
    :JumpProblem,
    # BipartiteGraphs
    :BipartiteGraph,
    # UnPack
    Symbol("@pack!"), Symbol("@unpack"), :UnPack,
)

run_qa(
    ModelingToolkitBase;
    Aqua,
    JET,
    jet = true,
    aqua_kwargs = (
        ;
        stale_deps = (; ignore = STALE_EXTENSION_DEPENDENCIES),
        piracies = (; treat_as_own = INTENTIONAL_EXTERNAL_GENERIC_EXTENSIONS),
    ),
    ei_kwargs = (;
        no_implicit_imports = (; allow_unanalyzable = GENERATED_ADT_MODULES),
        no_stale_explicit_imports = (; allow_unanalyzable = GENERATED_ADT_MODULES),
        all_explicit_imports_are_public = (; ignore = NONPUBLIC_EXPLICIT_IMPORTS),
        all_qualified_accesses_are_public = (; ignore = NONPUBLIC_QUALIFIED_ACCESSES),
    ),
    reexports_allow = REEXPORTED_API,
)
