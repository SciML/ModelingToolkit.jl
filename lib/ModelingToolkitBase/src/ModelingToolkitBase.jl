"""
$(DocStringExtensions.README)
"""
module ModelingToolkitBase

if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@max_methods"))
    @eval Base.Experimental.@compiler_options max_methods = 1
end

using PrecompileTools: PrecompileTools, @recompile_invalidations
using Reexport: Reexport, @reexport
@recompile_invalidations begin
    import StaticArrays
    using StaticArraysCore: StaticArraysCore, MVector, SVector, StaticArray, StaticVector,
        similar_type
    import Symbolics
    import ImplicitDiscreteSolve
    using ImplicitDiscreteSolve: IDSolve
    import JumpProcesses
    using JumpProcesses: ConstantRateJump, JumpProblem, JumpSet, MassActionJump,
        VariableRateJump, get_num_majumps, needs_depgraph, needs_vartojumps_map,
        reset_aggregated_jumps!
    # ONLY here for the invalidations
    import REPL
    using OffsetArrays: Origin
    import BlockArrays: BlockedArray, Block, blocksize, blocksizes
    import BandedMatrices: BandedMatrices, BandedMatrix, bandwidths
end

import SciMLBase
using SciMLBase: BVPFunction, BVProblem, CallbackSet, ContinuousCallback, DAEFunction,
    DAEProblem, DDEFunction, DDEProblem, DiscreteCallback, DiscreteFunction,
    DiscreteProblem, HomotopyNonlinearFunction, ImplicitDiscreteFunction,
    ImplicitDiscreteProblem, IntervalNonlinearFunction, IntervalNonlinearProblem,
    LinearProblem, NonlinearFunction, NonlinearLeastSquaresProblem, NonlinearProblem,
    ODEFunction, ODEInputFunction, ODEProblem, ODESolution, OptimizationFunction,
    OptimizationProblem, ReturnCode, SCCNonlinearProblem, SDDEFunction, SDDEProblem,
    SDEFunction, SDEProblem, SteadyStateProblem, VectorContinuousCallback, check_error,
    remake
using Printf: @sprintf

import SymbolicUtils
import SymbolicUtils as SU
import SymbolicUtils: iscall, arguments, operation, issym, FnType, isconst, BSImpl,
    @rule, Rewriters, substitute, BasicSymbolic, _iszero
using SymbolicUtils: @syms, BS, IRStructure, SymReal, expand, getmetadata, populate_ir!,
    setmetadata, simplify, simplify_fractions, unwrap_const
import TermInterface
import TermInterface: maketerm, metadata
import SymbolicUtils.Code
using SymbolicUtils.Code: Assignment, DestructuredArgs, Func, Let, MakeTuple, cse
import SymbolicUtils.Code: toexpr
import DocStringExtensions
using DocStringExtensions: FIELDS, METHODLIST, SIGNATURES, TYPEDEF, TYPEDFIELDS,
    TYPEDSIGNATURES
using SpecialFunctions: SpecialFunctions, gamma
import NaNMath
@recompile_invalidations begin
    import DiffEqCallbacks
    using DiffEqCallbacks: PeriodicCallback, PresetTimeCallback
    import DiffEqBase
    import ForwardDiff
end
import Graphs
using Graphs: BFSIterator, dfs_parents, dst, edges, simplecycles_iter, src,
    strongly_connected_components, topological_sort
import ExprTools: splitdef, combinedef
import OrderedCollections
using OrderedCollections: OrderedDict, OrderedSet

import SymbolicIndexingInterface
using SymbolicIndexingInterface: ArraySymbolic, ContinuousTimeseries, NotSymbolic,
    ParameterTimeseriesCollection, ParameterTimeseriesIndex, ProblemState, ScalarSymbolic,
    all_symbols, all_variable_symbols, current_time, getname, getsym, getu, hasname,
    independent_variable_symbols, is_independent_variable, is_markovian, is_parameter,
    is_time_dependent, is_timeseries_parameter, is_variable, parameter_index,
    parameter_symbols, parameter_values, set_parameter!, setp, setp_oop, setsym, setu,
    state_values, symbolic_container, symbolic_type, timeseries_parameter_index,
    variable_index, variable_symbols
import LinearAlgebra
using LinearAlgebra: Diagonal, I, UniformScaling, diag, diagm, dot, isdiag, tr
import SparseArrays
using SparseArrays: AbstractSparseArray, SparseMatrixCSC, findnz, nonzeros, sparse
import InteractiveUtils
import DataStructures
using DataStructures: Queue, dequeue!, enqueue!
using Base.Threads
import ArrayInterface
import Setfield
using Setfield: @set, @set!
import ConstructionBase
using ConstructionBase: constructorof, setproperties
import Libdl
import Combinatorics
import FunctionWrappersWrappers
import FunctionWrappersWrappers: FunctionWrappersWrapper
import FunctionWrappers: FunctionWrapper
import SciMLStructures
import Compat
import AbstractTrees
using AbstractTrees: TreeIterator, print_tree
using SciMLBase: StandardODEProblem, StandardNonlinearProblem, TimeDomain,
    Clock, SolverStepClock, AbstractNonlinearProblem
import Moshi
import RecursiveArrayTools
using RecursiveArrayTools: ArrayPartition, DiffEqArray
import Graphs: SimpleDiGraph, add_edge!
import CommonSolve
using CommonSolve: init, solve
import EnumX
import ReadOnlyDicts: ReadOnlyDict

import RuntimeGeneratedFunctions
using RuntimeGeneratedFunctions: RuntimeGeneratedFunction, drop_expr

using Symbolics: VartypeT, SymbolicT
using Symbolics: value, @derivatives, get_variables,
    symbolic_linear_solve, unwrap, wrap,
    VariableSource, variable,
    NAMESPACE_SEPARATOR, setdefaultval, Arr,
    fixpoint_sub, CallAndWrap, SArgsT, SSym, STerm, SConst
using Symbolics: @register_array_symbolic, @register_symbolic, @variables, Differential,
    Equation, Inequality, Integral, expand_derivatives, ≲, ≳
const NAMESPACE_SEPARATOR_SYMBOL = Symbol(NAMESPACE_SEPARATOR)
import Symbolics: rename, get_variables!,
    jacobian_sparsity, isaffine, islinear,
    tosymbol, lower_varname, diff2term, var_from_nested_derivative,
    build_function, linear_expansion, jacobian, sparsejacobian,
    scalarize, hasderiv

import SciMLBase: @add_kwonly
export independent_variables, unknowns, observables, parameters, bound_parameters,
    continuous_events, discrete_events, analytically_integrated
@reexport using Symbolics
@reexport using UnPack
import UnPack
using UnPack: @unpack
RuntimeGeneratedFunctions.init(@__MODULE__)

import DifferentiationInterface as DI
import SciMLPublic: @public
import PreallocationTools
import PreallocationTools: DiffCache, get_tmp
import FillArrays
import BipartiteGraphs
using BipartiteGraphs: BipartiteGraph, DiCMOBiGraph, HyperGraph, Matching, Unassigned,
    𝑑neighbors, 𝑠neighbors, 𝑠vertices
import Random: AbstractRNG
# To handle `Integral` in a type-stable manner
import DomainSets
# `DomainSets` re-exports IntervalSets' `endpoints`; reach it through its owner.
import IntervalSets
# For `LinearInitializationProblem`
import SCCNonlinearSolve
using TaskLocalValues: TaskLocalValue

export @derivatives

for fun in [:toexpr]
    @eval begin
        function $fun(eq::Equation; kw...)
            return Expr(:call, :(==), $fun(eq.lhs; kw...), $fun(eq.rhs; kw...))
        end

        function $fun(ineq::Inequality; kw...)
            return if ineq.relational_op == Symbolics.leq
                Expr(:call, :(<=), $fun(ineq.lhs; kw...), $fun(ineq.rhs; kw...))
            else
                Expr(:call, :(>=), $fun(ineq.lhs; kw...), $fun(ineq.rhs; kw...))
            end
        end

        $fun(eqs::AbstractArray; kw...) = map(eq -> $fun(eq; kw...), eqs)
        $fun(x::Integer; kw...) = x
        $fun(x::AbstractFloat; kw...) = x
    end
end

const INTERNAL_FIELD_WARNING = """
This field is internal API. It may be removed or changed without notice in a non-breaking \
release. Usage of this field is not advised.
"""

const INTERNAL_ARGS_WARNING = """
The following arguments are internal API. They may be removed or changed without notice \
in a non-breaking release. Usage of these arguments is not advised.
"""

"""
$(TYPEDEF)

Abstract supertype of all system types.

Custom system types must subtype `AbstractSystem` and implement the required structural
interface:

- `nameof(sys)::Symbol` must identify the system within its parent's subsystem list.
- `get_systems(sys)::Vector{<:AbstractSystem}` must return the direct, finite, acyclic
  subsystem hierarchy.

Both requirements default to fields named `name` and `systems`, respectively. Additional
system data is optional. Generic code must use the public `has_x`/`get_x` accessors described
in the [AbstractSystem interface](@ref abstract_system_interface), rather than depending on
the fields of `System` or another concrete subtype. A custom type may extend the documented
generic accessors, including `independent_variable`, when its storage does not match those
defaults; it should not extend internal compilation or problem-construction functions.

# Examples

A minimal custom system can participate in generic hierarchy traversal without implementing
any methods:

```julia
struct CustomSystem <: AbstractSystem
    name::Symbol
    systems::Vector{AbstractSystem}
end

leaf = CustomSystem(:leaf, AbstractSystem[])
root = CustomSystem(:root, AbstractSystem[leaf])
nameof(root)
get_systems(root)
```
"""
abstract type AbstractSystem end
# Solely so that `ODESystem` can be deprecated and still act as a valid type.
# See `deprecations.jl`.
abstract type IntermediateDeprecationSystem <: AbstractSystem end

"""
    independent_variable(sys)

Return the scalar independent variable of `sys`, or `nothing` when `sys` has no scalar
independent variable.

Most users should call [`independent_variables`](@ref), which always returns a vector. This
is the scalar extension point for packages that define custom `AbstractSystem` subtypes;
external packages should extend `independent_variable` only. The generic
`independent_variables` accessor wraps that scalar value and should not be extended.

# Arguments

- `sys`: A system-like object with a scalar independent variable.

# Returns

- The scalar independent variable of `sys`, or `nothing`.

# Examples

```julia
struct ScalarIVSystem <: ModelingToolkitBase.AbstractSystem
    iv
end
ModelingToolkitBase.independent_variable(sys::ScalarIVSystem) = getfield(sys, :iv)
```
"""
function independent_variable end
independent_variable(sys::AbstractSystem) = isdefined(sys, :iv) ? getfield(sys, :iv) : nothing

# this has to be included early to deal with dependency issues
function complete end

complete(m::Matching, args...; kw...) = BipartiteGraphs.complete(m, args...; kw...)
complete(g::BipartiteGraph, args...; kw...) = BipartiteGraphs.complete(g, args...; kw...)

export EvalAt
include("variables.jl")
include("parameters.jl")
include("discretes.jl")
include("independent_variables.jl")
include("constants.jl")
include("derivative_dict.jl")
include("atomic_array_dict.jl")
include("parameter_bindings_graph.jl")

const SymmapT = AtomicArrayDict{SymbolicT, Dict{SymbolicT, SymbolicT}}
const AtomicMapT{T} = AtomicArrayDict{T, Dict{SymbolicT, T}}
const AtomicSetT = AtomicArraySet{Dict{SymbolicT, Nothing}}
const ROSymmapT = ReadOnlyDict{SymbolicT, SymbolicT, SymmapT}
struct CommonSentinel end
const COMMON_SENTINEL = SU.Const{VartypeT}(CommonSentinel())
const COMMON_NOTHING = SU.Const{VartypeT}(nothing)
const COMMON_MISSING = SU.Const{VartypeT}(missing)
const COMMON_TRUE = SU.Const{VartypeT}(true)
const COMMON_FALSE = SU.Const{VartypeT}(false)
const COMMON_INF = SU.Const{VartypeT}(Inf)

include("utils.jl")

include("systems/index_cache.jl")
include("systems/parameter_buffer.jl")
include("systems/abstractsystem.jl")
# codegen_utils.jl defines the codegen option structs (GeneratedFunctionOptions,
# BuildFunctionWrapperOptions, CompilerOptions). It must precede any file whose method
# *signatures* annotate those types (e.g. the `opts::GeneratedFunctionOptions` positional
# in callbacks.jl), since signature type annotations are resolved at definition time.
include("systems/codegen_utils.jl")
include("systems/connectiongraph.jl")
include("systems/connectors.jl")
include("systems/imperative_affect.jl")
include("systems/callbacks.jl")
include("systems/system.jl")
include("systems/analysis_points.jl")
include("systems/ir_info.jl")
include("problems/docs.jl")
include("systems/codegen.jl")
include("systems/codegen_compat.jl")
include("systems/problem_utils.jl")
# Operator + lowering layer; must load before the problem constructors that
# consume it (problems/nonlinearproblem.jl selector + problems/homotopyproblem.jl).
include("systems/homotopy_operator.jl")

include("problems/compatibility.jl")
include("problems/odeproblem.jl")
include("problems/ddeproblem.jl")
include("problems/daeproblem.jl")
include("problems/sdeproblem.jl")
include("problems/sddeproblem.jl")
include("problems/nonlinearproblem.jl")
include("problems/homotopyproblem.jl")
include("problems/intervalnonlinearproblem.jl")
include("problems/implicitdiscreteproblem.jl")
include("problems/discreteproblem.jl")
include("problems/optimizationproblem.jl")
include("problems/jumpproblem.jl")
include("problems/initializationproblem.jl")
include("problems/bvproblem.jl")
include("problems/linearproblem.jl")

include("modelingtoolkitize/common.jl")
include("modelingtoolkitize/odeproblem.jl")
include("modelingtoolkitize/sdeproblem.jl")
include("modelingtoolkitize/optimizationproblem.jl")
include("modelingtoolkitize/nonlinearproblem.jl")

include("systems/nonlinear/homotopy_continuation.jl")
include("systems/nonlinear/initializesystem.jl")
include("systems/diffeqs/basic_transformations.jl")

include("systems/pde/pdesystem.jl")


include("systems/unit_check.jl")
include("systems/dependency_graphs.jl")
include("discretedomain.jl")
include("systems/systems.jl")

include("debugging.jl")

include("inputoutput.jl")

include("deprecations.jl")

"""
    t_nounits

Unitless default independent variable used by ModelingToolkit examples and constructors.

# Examples

```julia
using ModelingToolkitBase

t = ModelingToolkitBase.t_nounits
```
"""
const t_nounits = let
    only(@independent_variables t)
end

"""
    D_nounits

Default unitless differential operator `Differential(t_nounits)`.

# Examples

```julia
using ModelingToolkitBase

D = ModelingToolkitBase.D_nounits
```
"""
const D_nounits = Differential(t_nounits)

export CompilerOptions
export ODEFunction, convert_system_indepvar,
    System, OptimizationSystem, JumpSystem, SDESystem, NonlinearSystem, ODESystem,
    DiscreteSystem, ImplicitDiscreteSystem
export SDEFunction
export DiscreteProblem, DiscreteFunction
export ImplicitDiscreteProblem, ImplicitDiscreteFunction
export ODEProblem, SDEProblem
export NonlinearFunction
export NonlinearProblem
# `AbstractNonlinearProblem(sys, op)` is the automatic constructor that selects a
# `HomotopyProblem` when `sys` contains `homotopy(...)` nodes (see
# `problems/homotopyproblem.jl`); re-exported so it is reachable unqualified.
export AbstractNonlinearProblem
export IntervalNonlinearFunction
export IntervalNonlinearProblem
export OptimizationProblem, constraints
export SteadyStateProblem
export JumpProblem, SymbolicMassActionJump
export flatten
export connect, domain_connect, @connector, Connection, AnalysisPoint, Flow, Stream,
    instream
export @component, @mtkcompile, @mtkbuild, @mtkcomplete
export isinput, isoutput, getbounds, hasbounds, getnominal, hasnominal, setnominal, getguess, hasguess, isdisturbance,
    istunable, getdist, hasdist,
    tunable_parameters, isirreducible, getdescription, hasdescription,
    hasunit, getunit, hasconnect, getconnect,
    hasmisc, getmisc, state_priority,
    subset_tunables
export liouville_transform, change_independent_variable,
    add_accumulations, noise_to_brownians, Girsanov_transform, change_of_variables,
    fractional_to_ordinary, linear_fractional_to_ordinary
export respecialize
export PDESystem
export Differential, expand_derivatives, @derivatives
export Equation
export Term
export SymScope, LocalScope, ParentScope, GlobalScope
export independent_variable, equations, observed, full_equations, jumps, cost,
    brownians
export initialization_equations, guesses, bindings, initial_conditions, hierarchy
export set_defaults
export state_priorities, irreducibles, maybe_zeros
export mtkcompile, expand_connections, structural_simplify
export solve
export Pre

export calculate_jacobian, generate_jacobian, generate_rhs, generate_custom_function,
    generate_W, calculate_hessian
export calculate_control_jacobian, generate_control_jacobian
export calculate_tgrad, generate_tgrad
export generate_cost, calculate_cost_gradient, generate_cost_gradient
export generate_trajectory
export calculate_cost_hessian, generate_cost_hessian
export calculate_massmatrix, generate_diffusion_function
export generate_control_function, build_explicit_observed_function
export stochastic_integral_transform

export BipartiteGraph, equation_dependencies, variable_dependencies
export eqeq_dependencies, varvar_dependencies
export asgraph, asdigraph

export toexpr, get_variables
export simplify, substitute
export build_function
export modelingtoolkitize
export generate_initializesystem, Initial, isinitial, InitializationProblem

export alg_equations, diff_equations, has_alg_equations, has_diff_equations
export get_alg_eqs, get_diff_eqs, has_alg_eqs, has_diff_eqs

export homotopy

export @variables, @parameters, @independent_variables, @constants, @brownians, @brownian,
    @poissonians, @discretes
export @named, @nonamespace, @namespace, extend, compose, complete, toggle_namespacing
export debug_system

#export ContinuousClock, Discrete, sampletime, input_timedomain, output_timedomain
#export has_discrete_domain, has_continuous_domain
#export is_discrete_domain, is_continuous_domain, is_hybrid_domain
export Shift, ShiftIndex
export Sample, Hold, SampleTime
export Clock, SolverStepClock, TimeDomain

export MTKParameters, reorder_dimension_by_tunables!, reorder_dimension_by_tunables

export HomotopyContinuationProblem

export AnalysisPoint, open_loop

include("systems/optimal_control_interface.jl")

using SciMLBase: AbstractDynamicOptProblem
export AbstractDynamicOptProblem, JuMPDynamicOptProblem, InfiniteOptDynamicOptProblem,
    CasADiDynamicOptProblem, PyomoDynamicOptProblem
export AbstractCollocation, JuMPCollocation, InfiniteOptCollocation,
    CasADiCollocation, PyomoCollocation
export DynamicOptSolution
export MissingGuessValue
export MiscSystemData

export AssignmentAffect

const set_scalar_metadata = setmetadata

@public apply_to_variables, equations_toplevel, unknowns_toplevel, parameters_toplevel
@public continuous_events_toplevel, discrete_events_toplevel, assertions, is_alg_equation
@public is_diff_equation, Equality
@public inputs, outputs, bound_inputs, unbound_inputs, bound_outputs
@public unbound_outputs, is_bound
@public AbstractSystem, CheckAll, CheckNone, CheckComponents, CheckUnits
@public t, D, t_nounits, D_nounits
@public SymbolicContinuousCallback, SymbolicDiscreteCallback, ImperativeAffect
@public VariableType, MTKVariableTypeCtx, VariableBounds, VariableConnectType
@public VariableDescription, VariableInput, VariableIrreducible, VariableMisc
@public VariableOutput, VariableStatePriority, VariableUnit, collect_scoped_vars!
@public collect_var_to_name!, collect_vars!, eqtype_supports_collect_vars, hasdefault
@public getdefault, setdefault, setguess, iscomplete, isparameter, modified_unknowns!
@public renamespace, namespace_equations
@public varmap_to_vars
@public check_mutable_cache, store_to_mutable_cache!, should_invalidate_mutable_cache_entry
@public convert_bindings_for_time_independent_system, get_w
@public Both
@public SymbolicADDisallowed, check_symbolic_ad_allowed
@public tobrownian, toparam
@public ProblemTypeCtx

for prop in [SYS_PROPS; [:continuous_events, :discrete_events]]
    getter = Symbol(:get_, prop)
    hasfn = Symbol(:has_, prop)
    @eval @public $getter, $hasfn
end

# AbstractSystem types should be treated as inactive (constant) for Enzyme.
# Property access on systems retrieves symbolic metadata, not numerical values.
import EnzymeCore
function EnzymeCore.EnzymeRules.inactive_noinl(
        ::typeof(Base.getproperty), ::AbstractSystem, ::Symbol,
    )
    return true
end
# Declare the type itself inactive so Enzyme's runtime-activity dispatch
# (e.g. `MixedDuplicated`-wrapping) never tries to track a `System` as a
# differentiable value. Without this, `Enzyme.gradient(set_runtime_activity(Reverse),
# Const(loss), p)` over a closure capturing an MTK-generated problem
# (which transitively carries the symbolic `System`) trips a
# `MethodError MixedDuplicated(::System, ::System)` in `create_activity_wrapper`.
EnzymeCore.EnzymeRules.inactive_type(::Type{<:AbstractSystem}) = true

function __init__()
    SU.hashcons(unwrap(t_nounits), true)
    SU.hashcons(COMMON_NOTHING, true)
    SU.hashcons(COMMON_MISSING, true)
    SU.hashcons(COMMON_TRUE, true)
    SU.hashcons(COMMON_FALSE, true)
    SU.hashcons(COMMON_SENTINEL, true)
    SU.hashcons(COMMON_INF, true)
    SU.hashcons(MTKPARAMETERS_ARG, true)
    SU.hashcons(COMMON_DEFAULT_VAR, true)
    SU.hashcons(DDE_HISTORY_FUN, true)
    SU.hashcons(BVP_SOLUTION, true)
    SU.hashcons(unwrap(W_GAMMA), true)
    SU.hashcons(unwrap(ASSERTION_LOG_VARIABLE), true)
    SU.hashcons(DDE_AT_IDX_SYM, true)
    SU.hashcons(DDE_DELAY_SYM, true)
    SU.hashcons(MTKUNKNOWNS_ARG, true)
    return nothing
end

include("precompile.jl")
end # module
