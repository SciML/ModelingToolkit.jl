"""
$(DocStringExtensions.README)
"""
module ModelingToolkit
using PrecompileTools, Reexport
@recompile_invalidations begin
    using StaticArrays
    using Symbolics
end

import SymbolicUtils
import SymbolicUtils: iscall, arguments, operation, maketerm, promote_symtype,
                      Symbolic, isadd, ismul, ispow, issym, FnType,
                      @rule, Rewriters, substitute, metadata, BasicSymbolic,
                      Sym, Term
using SymbolicUtils.Code
import SymbolicUtils.Code: toexpr
import SymbolicUtils.Rewriters: Chain, Postwalk, Prewalk, Fixpoint
using DocStringExtensions
using SpecialFunctions, NaNMath
using DiffEqCallbacks
using Graphs
import ExprTools: splitdef, combinedef
import OrderedCollections
using DiffEqNoiseProcess: DiffEqNoiseProcess, WienerProcess

using SymbolicIndexingInterface
using LinearAlgebra, SparseArrays
using InteractiveUtils
using JumpProcesses
using DataStructures
@static if pkgversion(DataStructures) >= v"0.19"
    import DataStructures: IntDisjointSet
else
    import DataStructures: IntDisjointSets
    const IntDisjointSet = IntDisjointSets
end
using Base.Threads
using Latexify, Unitful, ArrayInterface
using Setfield, ConstructionBase
import Libdl
using DocStringExtensions
using Base: RefValue
using Combinatorics
import Distributions
import FunctionWrappersWrappers
import FunctionWrappers: FunctionWrapper
using URIs: URI
using SciMLStructures
using Compat
using AbstractTrees
using DiffEqBase, SciMLBase, ForwardDiff
using SciMLBase: StandardODEProblem, StandardNonlinearProblem, handle_varmap, TimeDomain,
                 PeriodicClock, Clock, SolverStepClock, ContinuousClock, OverrideInit,
                 NoInit
using Distributed
import JuliaFormatter
using MLStyle
import Moshi
using Moshi.Data: @data
import SCCNonlinearSolve
using ImplicitDiscreteSolve
using Reexport
using RecursiveArrayTools
import Graphs: SimpleDiGraph, add_edge!, incidence_matrix
import BlockArrays: BlockArray, BlockedArray, Block, blocksize, blocksizes, blockpush!,
                    undef_blocks, blocks
using OffsetArrays: Origin
import CommonSolve
import EnumX
import ChainRulesCore
import ChainRulesCore: Tangent, ZeroTangent, NoTangent, zero_tangent, unthunk

using RuntimeGeneratedFunctions
using RuntimeGeneratedFunctions: drop_expr

using Symbolics: degree
using Symbolics: _parse_vars, value, @derivatives, get_variables,
                 exprs_occur_in, symbolic_linear_solve, build_expr, unwrap, wrap,
                 VariableSource, getname, variable,
                 NAMESPACE_SEPARATOR, set_scalar_metadata, setdefaultval,
                 hasnode, fixpoint_sub, fast_substitute,
                 CallWithMetadata, CallWithParent
const NAMESPACE_SEPARATOR_SYMBOL = Symbol(NAMESPACE_SEPARATOR)
import Symbolics: rename, get_variables!, _solve, hessian_sparsity,
                  jacobian_sparsity, isaffine, islinear, _iszero, _isone,
                  tosymbol, lower_varname, diff2term, var_from_nested_derivative,
                  BuildTargets, JuliaTarget, StanTarget, CTarget, MATLABTarget,
                  ParallelForm, SerialForm, MultithreadedForm, build_function,
                  rhss, lhss, prettify_expr, gradient,
                  jacobian, hessian, derivative, sparsejacobian, sparsehessian,
                  substituter, scalarize, getparent, hasderiv, hasdiff

import DiffEqBase: @add_kwonly
export independent_variables, unknowns, observables, parameters, full_parameters,
       continuous_events, discrete_events
@reexport using Symbolics
@reexport using UnPack
RuntimeGeneratedFunctions.init(@__MODULE__)

import DynamicQuantities, Unitful
const DQ = DynamicQuantities

import DifferentiationInterface as DI
using ADTypes: AutoForwardDiff
import SciMLPublic: @public

export @derivatives

for fun in [:toexpr]
    @eval begin
        function $fun(eq::Equation; kw...)
            Expr(:call, :(==), $fun(eq.lhs; kw...), $fun(eq.rhs; kw...))
        end

        function $fun(ineq::Inequality; kw...)
            if ineq.relational_op == Symbolics.leq
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

function independent_variable end

# this has to be included early to deal with dependency issues
include("structural_transformation/bareiss.jl")
function complete end
function var_derivative! end
function var_derivative_graph! end
include("bipartite_graph.jl")
using .BipartiteGraphs

export EvalAt
include("variables.jl")
include("parameters.jl")
include("independent_variables.jl")
include("constants.jl")

include("utils.jl")

include("systems/index_cache.jl")
include("systems/parameter_buffer.jl")
include("systems/abstractsystem.jl")
include("systems/model_parsing.jl")
include("systems/connectiongraph.jl")
include("systems/connectors.jl")
include("systems/state_machines.jl")
include("systems/analysis_points.jl")
include("systems/imperative_affect.jl")
include("systems/callbacks.jl")
include("systems/system.jl")
include("systems/codegen_utils.jl")
include("problems/docs.jl")
include("systems/codegen.jl")
include("systems/problem_utils.jl")
include("linearization.jl")
include("systems/solver_nlprob.jl")

include("problems/compatibility.jl")
include("problems/odeproblem.jl")
include("problems/ddeproblem.jl")
include("problems/daeproblem.jl")
include("problems/sdeproblem.jl")
include("problems/sddeproblem.jl")
include("problems/nonlinearproblem.jl")
include("problems/intervalnonlinearproblem.jl")
include("problems/implicitdiscreteproblem.jl")
include("problems/discreteproblem.jl")
include("problems/optimizationproblem.jl")
include("problems/jumpproblem.jl")
include("problems/initializationproblem.jl")
include("problems/sccnonlinearproblem.jl")
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

include("systems/sparsematrixclil.jl")

include("systems/unit_check.jl")
include("systems/validation.jl")
include("systems/dependency_graphs.jl")
include("clock.jl")
include("discretedomain.jl")
include("systems/systemstructure.jl")
include("systems/clock_inference.jl")
include("systems/systems.jl")
include("systems/if_lifting.jl")

include("debugging.jl")
include("systems/alias_elimination.jl")
include("structural_transformation/StructuralTransformations.jl")

@reexport using .StructuralTransformations
include("inputoutput.jl")

include("adjoints.jl")
include("deprecations.jl")

const t_nounits = let
    only(@independent_variables t)
end
const t_unitful = let
    only(@independent_variables t [unit = Unitful.u"s"])
end
const t = let
    only(@independent_variables t [unit = DQ.u"s"])
end

const D_nounits = Differential(t_nounits)
const D_unitful = Differential(t_unitful)
const D = Differential(t)

export ODEFunction, convert_system_indepvar,
       System, OptimizationSystem, JumpSystem, SDESystem, NonlinearSystem, ODESystem
export SDEFunction
export SystemStructure
export DiscreteProblem, DiscreteFunction
export ImplicitDiscreteProblem, ImplicitDiscreteFunction
export ODEProblem, SDEProblem
export NonlinearFunction
export NonlinearProblem
export IntervalNonlinearFunction
export IntervalNonlinearProblem
export OptimizationProblem, constraints
export SteadyStateProblem
export JumpProblem
export alias_elimination, flatten
export connect, domain_connect, @connector, Connection, AnalysisPoint, Flow, Stream,
       instream
export initial_state, transition, activeState, entry, ticksInState, timeInState
export @component, @mtkmodel, @mtkcompile, @mtkbuild
export isinput, isoutput, getbounds, hasbounds, getguess, hasguess, isdisturbance,
       istunable, getdist, hasdist,
       tunable_parameters, isirreducible, getdescription, hasdescription,
       hasunit, getunit, hasconnect, getconnect,
       hasmisc, getmisc, state_priority,
       subset_tunables
export liouville_transform, change_independent_variable, substitute_component,
       add_accumulations, noise_to_brownians, Girsanov_transform, change_of_variables,
       fractional_to_ordinary, linear_fractional_to_ordinary
export respecialize
export PDESystem
export Differential, expand_derivatives, @derivatives
export Equation, ConstrainedEquation
export Term, Sym
export SymScope, LocalScope, ParentScope, GlobalScope
export independent_variable, equations, observed, full_equations, jumps, cost,
       brownians
export initialization_equations, guesses, defaults, parameter_dependencies, hierarchy
export mtkcompile, expand_connections, linearize, linearization_function,
       LinearizationProblem, linearization_ap_transform, structural_simplify
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
export TearingState

export BipartiteGraph, equation_dependencies, variable_dependencies
export eqeq_dependencies, varvar_dependencies
export asgraph, asdigraph
export map_variables_to_equations

export toexpr, get_variables
export simplify, substitute
export build_function
export modelingtoolkitize
export generate_initializesystem, Initial, isinitial, InitializationProblem

export alg_equations, diff_equations, has_alg_equations, has_diff_equations
export get_alg_eqs, get_diff_eqs, has_alg_eqs, has_diff_eqs

export @variables, @parameters, @independent_variables, @constants, @brownians, @brownian
export @named, @nonamespace, @namespace, extend, compose, complete, toggle_namespacing
export debug_system

#export ContinuousClock, Discrete, sampletime, input_timedomain, output_timedomain
#export has_discrete_domain, has_continuous_domain
#export is_discrete_domain, is_continuous_domain, is_hybrid_domain
export Sample, Hold, Shift, ShiftIndex, sampletime, SampleTime
export Clock, SolverStepClock, TimeDomain

export MTKParameters, reorder_dimension_by_tunables!, reorder_dimension_by_tunables

export HomotopyContinuationProblem

export AnalysisPoint, get_sensitivity_function, get_comp_sensitivity_function,
       get_looptransfer_function, get_sensitivity, get_comp_sensitivity, get_looptransfer,
       open_loop
function FMIComponent end

include("systems/optimal_control_interface.jl")
export AbstractDynamicOptProblem, JuMPDynamicOptProblem, InfiniteOptDynamicOptProblem,
       CasADiDynamicOptProblem, PyomoDynamicOptProblem
export AbstractCollocation, JuMPCollocation, InfiniteOptCollocation,
       CasADiCollocation, PyomoCollocation
export DynamicOptSolution

@public apply_to_variables, equations_toplevel, unknowns_toplevel, parameters_toplevel
@public continuous_events_toplevel, discrete_events_toplevel, assertions, is_alg_equation
@public is_diff_equation, Equality, linearize_symbolic, reorder_unknowns
@public similarity_transform, inputs, outputs, bound_inputs, unbound_inputs, bound_outputs
@public unbound_outputs, is_bound
@public AbstractSystem, CheckAll, CheckNone, CheckComponents, CheckUnits
@public t, D, t_nounits, D_nounits, t_unitful, D_unitful
@public SymbolicContinuousCallback, SymbolicDiscreteCallback
@public VariableType, MTKVariableTypeCtx, VariableBounds, VariableConnectType
@public VariableDescription, VariableInput, VariableIrreducible, VariableMisc
@public VariableOutput, VariableStatePriority, VariableUnit, collect_scoped_vars!
@public collect_var_to_name!, collect_vars!, eqtype_supports_collect_vars, hasdefault
@public getdefault, setdefault, iscomplete, isparameter, modified_unknowns!
@public renamespace, namespace_equations

for prop in [SYS_PROPS; [:continuous_events, :discrete_events]]
    getter = Symbol(:get_, prop)
    hasfn = Symbol(:has_, prop)
    @eval @public $getter, $hasfn
end

PrecompileTools.@compile_workload begin
    using ModelingToolkit
    @variables x(ModelingToolkit.t_nounits)
    @named sys = System([ModelingToolkit.D_nounits(x) ~ -x], ModelingToolkit.t_nounits)
    prob = ODEProblem(mtkcompile(sys), [x => 30.0], (0, 100), jac = true)
    @mtkmodel __testmod__ begin
        @constants begin
            c = 1.0
        end
        @structural_parameters begin
            structp = false
        end
        if structp
            @variables begin
                x(t) = 0.0, [description = "foo", guess = 1.0]
            end
        else
            @variables begin
                x(t) = 0.0, [description = "foo w/o structp", guess = 1.0]
            end
        end
        @parameters begin
            a = 1.0, [description = "bar"]
            if structp
                b = 2 * a, [description = "if"]
            else
                c
            end
        end
        @equations begin
            x ~ a + b
        end
    end
end

end # module
