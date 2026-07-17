"""
$(DocStringExtensions.README)
"""
module ModelingToolkitBase

if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@max_methods"))
    @eval Base.Experimental.@compiler_options max_methods = 1
end

using PrecompileTools, Reexport
@recompile_invalidations begin
    using StaticArrays
    using Symbolics
    using ImplicitDiscreteSolve
    using JumpProcesses
    # ONLY here for the invalidations
    import REPL
    using OffsetArrays: Origin
    import BlockArrays: BlockArray, BlockedArray, Block, blocksize, blocksizes, blockpush!,
        undef_blocks, blocks
    import BandedMatrices: BandedMatrices, BandedMatrix, bandwidths
end

import SymbolicUtils
import SymbolicUtils as SU
import SymbolicUtils: iscall, arguments, operation, maketerm, promote_symtype,
    isadd, ismul, ispow, issym, FnType, isconst, BSImpl,
    @rule, Rewriters, substitute, metadata, BasicSymbolic
import SymbolicUtils: Code
using SymbolicUtils.Code: Assignment, DestructuredArgs, Func, Let, MakeTuple
import SymbolicUtils.Rewriters: Chain, Postwalk, Prewalk, Fixpoint
using DocStringExtensions
using SpecialFunctions, NaNMath
@recompile_invalidations begin
    using DiffEqCallbacks
    using DiffEqBase, SciMLBase, ForwardDiff
end
using Graphs
import ExprTools: splitdef, combinedef
import OrderedCollections

using SymbolicIndexingInterface
using LinearAlgebra, SparseArrays
using InteractiveUtils
using DataStructures
using Base.Threads
using ArrayInterface
using Setfield, ConstructionBase
import Libdl
using DocStringExtensions
using Base: RefValue
using Combinatorics
import FunctionWrappersWrappers
import FunctionWrappersWrappers: FunctionWrappersWrapper
import FunctionWrappers: FunctionWrapper
using SciMLStructures
using Compat
using AbstractTrees
using SciMLBase: StandardODEProblem, StandardNonlinearProblem, handle_varmap, TimeDomain,
    PeriodicClock, Clock, SolverStepClock, ContinuousClock, OverrideInit,
    NoInit, AbstractNonlinearProblem
import Moshi
using Moshi.Data: @data
using Reexport
using RecursiveArrayTools
import Graphs: SimpleDiGraph, add_edge!, incidence_matrix
import CommonSolve
import EnumX
import ReadOnlyDicts: ReadOnlyDict

using RuntimeGeneratedFunctions
using RuntimeGeneratedFunctions: drop_expr

using Symbolics: degree, VartypeT, SymbolicT
using Symbolics: parse_vars, value, @derivatives, get_variables,
    exprs_occur_in, symbolic_linear_solve, unwrap, wrap,
    VariableSource, getname, variable,
    NAMESPACE_SEPARATOR, setdefaultval, Arr,
    hasnode, fixpoint_sub, CallAndWrap, SArgsT, SSym, STerm, SConst
const NAMESPACE_SEPARATOR_SYMBOL = Symbol(NAMESPACE_SEPARATOR)
import Symbolics: rename, get_variables!, _solve, hessian_sparsity,
    jacobian_sparsity, isaffine, islinear, _iszero, _isone,
    tosymbol, lower_varname, diff2term, var_from_nested_derivative,
    BuildTargets, JuliaTarget, StanTarget, CTarget, MATLABTarget,
    ParallelForm, SerialForm, MultithreadedForm, build_function,
    rhss, lhss, gradient, linear_expansion,
    jacobian, hessian, derivative, sparsejacobian, sparsehessian,
    scalarize, hasderiv

import DiffEqBase: @add_kwonly
export independent_variables, unknowns, observables, parameters, bound_parameters,
    continuous_events, discrete_events, analytically_integrated
@reexport using Symbolics
@reexport using UnPack
RuntimeGeneratedFunctions.init(@__MODULE__)

"""
    @symbolic_wrap expr

Evaluate `expr` and wrap any symbolic result in the high-level Symbolics wrapper type.

This is a Symbolics utility reexported by ModelingToolkitBase. It is mainly useful when
writing low-level symbolic utilities that must preserve the wrapped `Num` interface.

# Examples

```julia
using ModelingToolkitBase

@variables x
wrapped_x = @symbolic_wrap Symbolics.unwrap(x)
```
"""
var"@symbolic_wrap"

"""
    @wrapped expr
    @wrapped expr wrap_arrays

Evaluate `expr` with symbolic arguments unwrapped, then wrap the result.

This is a Symbolics utility reexported by ModelingToolkitBase. It is useful for defining
helpers that operate on the internal symbolic representation while returning wrapped
values to users.

# Examples

```julia
using ModelingToolkitBase

@variables x
y = @wrapped Symbolics.unwrap(x)
```
"""
var"@wrapped"

"""
    RuleSet(rules)

Container for a collection of SymbolicUtils rewrite rules.

`RuleSet` is reexported for users who build custom symbolic simplification or rewriting
passes before handing expressions to ModelingToolkit.

# Examples

```julia
using ModelingToolkitBase
using SymbolicUtils: @rule

rules = RuleSet([@rule sin(~x)^2 + cos(~x)^2 => 1])
```
"""
RuleSet

"""
    get_canonical_expr(ir, idx)

Return the canonical expression stored at index `idx` in a SymbolicUtils IR structure.

This is a low-level SymbolicUtils helper reexported through Symbolics. Most ModelingToolkit
users should work with symbolic expressions directly instead of using the IR interface.

# Examples

```julia
using ModelingToolkitBase

# `get_canonical_expr(ir, idx)` is used after constructing a SymbolicUtils IR.
```
"""
get_canonical_expr

"""
    infimum(domain)

Return the lower endpoint of a symbolic interval domain.

This extends DomainSets interval queries to intervals whose endpoints are Symbolics
expressions.

# Examples

```julia
using DomainSets, ModelingToolkitBase

@variables a b
lo = infimum(a .. b)
```
"""
infimum

"""
    supremum(domain)

Return the upper endpoint of a symbolic interval domain.

This extends DomainSets interval queries to intervals whose endpoints are Symbolics
expressions.

# Examples

```julia
using DomainSets, ModelingToolkitBase

@variables a b
hi = supremum(a .. b)
```
"""
supremum

"""
    is_derivative(x) -> Bool

Return whether `x` is a symbolic derivative expression.

# Examples

```julia
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D

@variables x(t)
is_derivative(Symbolics.unwrap(D(x)))
```
"""
is_derivative

"""
    istree(x) -> Bool

Deprecated alias for `iscall(x)` from SymbolicUtils.

Use `iscall` in new code to check whether a symbolic value has an operation and
arguments.

# Examples

```julia
using ModelingToolkitBase

@variables x
iscall(Symbolics.unwrap(sin(x)))
```
"""
istree

"""
    solve_for(eq, var; simplify = false, check = true)

Solve a symbolic equation or equation system for `var`.

# Examples

```julia
using ModelingToolkitBase

@variables x y
solve_for(2x + y ~ 0, x)
```
"""
solve_for

import DifferentiationInterface as DI
using ADTypes: AutoForwardDiff
import SciMLPublic: @public
import PreallocationTools
import PreallocationTools: DiffCache, get_tmp
import FillArrays
using BipartiteGraphs
import Random: AbstractRNG
# To handle `Integral` in a type-stable manner
import DomainSets
# For `LinearInitializationProblem`
import SCCNonlinearSolve
using TaskLocalValues: TaskLocalValue

export @derivatives

"""
    toexpr(x; kwargs...)

Convert a symbolic ModelingToolkit object to a Julia expression.

`toexpr` is used by ModelingToolkit's generated-code and serialization utilities. For
plain SymbolicUtils code-generation objects, this delegates to `SymbolicUtils.Code.toexpr`.

# Examples

```julia
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t

@variables x(t)
expr = toexpr(Differential(t)(x) ~ x)
```
"""
toexpr(x; kw...) = Code.toexpr(x; kw...)

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

Abstract supertype of all system types. Any custom system types must subtype this.
"""
abstract type AbstractSystem end
# Solely so that `ODESystem` can be deprecated and still act as a valid type.
# See `deprecations.jl`.
abstract type IntermediateDeprecationSystem <: AbstractSystem end

"""
    independent_variable

Generic function for querying the primary independent variable of a system-like object.

Most users should call [`independent_variables`](@ref), which returns the independent
variables as a vector. Packages that define custom system types may extend
`independent_variable` when a scalar independent-variable interface is required.
"""
function independent_variable end

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
    System, OptimizationSystem, JumpSystem, SDESystem, NonlinearSystem, ODESystem
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
export state_priorities, irreducibles, maybe_zeros
export mtkcompile, expand_connections, structural_simplify
export solve
export Pre

export calculate_jacobian, generate_jacobian, generate_rhs, generate_custom_function,
    generate_W, calculate_hessian
export calculate_control_jacobian, generate_control_jacobian
export calculate_tgrad, generate_tgrad
export generate_cost, calculate_cost_gradient, generate_cost_gradient
export calculate_cost_hessian, generate_cost_hessian
export calculate_massmatrix, generate_diffusion_function
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
@public check_mutable_cache, store_to_mutable_cache!, should_invalidate_mutable_cache_entry
@public convert_bindings_for_time_independent_system, get_w
@public Both
@public SymbolicADDisallowed, check_symbolic_ad_allowed

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
