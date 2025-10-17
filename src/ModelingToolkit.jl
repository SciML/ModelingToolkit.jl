""""""
module ModelingToolkit
using PrecompileTools, Reexport
@recompile_invalidations begin
    using StaticArrays
    using Symbolics
    using ImplicitDiscreteSolve
    using JumpProcesses
    import REPL
end
import SymbolicUtils
import SymbolicUtils as SU
import SymbolicUtils: iscall, arguments, operation, maketerm, promote_symtype,
                      isadd, ismul, ispow, issym, FnType, isconst, BSImpl,
                      @rule, Rewriters, substitute, metadata, BasicSymbolic
using SymbolicUtils.Code
import SymbolicUtils.Code: toexpr
import SymbolicUtils.Rewriters: Chain, Postwalk, Prewalk, Fixpoint
using DocStringExtensions
using SpecialFunctions, NaNMath
@recompile_invalidations begin
    using DiffEqCallbacks
using DiffEqNoiseProcess: DiffEqNoiseProcess, WienerProcess
using DiffEqBase, SciMLBase, ForwardDiff
end
using Graphs
import ExprTools: splitdef, combinedef
import OrderedCollections
using SymbolicIndexingInterface
using LinearAlgebra, SparseArrays
using InteractiveUtils
using DataStructures
@static if pkgversion(DataStructures) >= v"0.19"
    import DataStructures: IntDisjointSet
else
    import DataStructures: IntDisjointSets
    const IntDisjointSet = IntDisjointSets
end
using Base.Threads
using Latexify, ArrayInterface
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
using SciMLBase: StandardODEProblem, StandardNonlinearProblem, handle_varmap, TimeDomain,
                 PeriodicClock, Clock, SolverStepClock, ContinuousClock, OverrideInit,
                 NoInit
using Distributed
using MLStyle
import Moshi
using Moshi.Data: @data
import SCCNonlinearSolve
using Reexport
using RecursiveArrayTools
import Graphs: SimpleDiGraph, add_edge!, incidence_matrix
import BlockArrays: BlockArray, BlockedArray, Block, blocksize, blocksizes, blockpush!,
                    undef_blocks, blocks
using OffsetArrays: Origin
import CommonSolve
import EnumX
using RuntimeGeneratedFunctions
using RuntimeGeneratedFunctions: drop_expr
using Symbolics: degree, VartypeT, SymbolicT
using Symbolics: parse_vars, value, @derivatives, get_variables,
                 exprs_occur_in, symbolic_linear_solve, unwrap, wrap,
                 VariableSource, getname, variable,
                 NAMESPACE_SEPARATOR, setdefaultval, Arr,
                 hasnode, fixpoint_sub, CallAndWrap, SArgsT, SSym, STerm
const NAMESPACE_SEPARATOR_SYMBOL = Symbol(NAMESPACE_SEPARATOR)
import Symbolics: rename, get_variables!, _solve, hessian_sparsity,
                  jacobian_sparsity, isaffine, islinear, _iszero, _isone,
                  tosymbol, lower_varname, diff2term, var_from_nested_derivative,
                  BuildTargets, JuliaTarget, StanTarget, CTarget, MATLABTarget,
                  ParallelForm, SerialForm, MultithreadedForm, build_function,
                  rhss, lhss, prettify_expr, gradient,
                  jacobian, hessian, derivative, sparsejacobian, sparsehessian,
                  scalarize, hasderiv
import DiffEqBase: @add_kwonly
export independent_variables, unknowns, observables, parameters, full_parameters,
       continuous_events, discrete_events
@reexport using Symbolics
@reexport using UnPack
RuntimeGeneratedFunctions.init(@__MODULE__)
import DynamicQuantities
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
""""""
abstract type AbstractSystem end
abstract type IntermediateDeprecationSystem <: AbstractSystem end
function independent_variable end
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
const SymmapT = Dict{SymbolicT, SymbolicT}
const COMMON_NOTHING = SU.Const{VartypeT}(nothing)
const COMMON_MISSING = SU.Const{VartypeT}(missing)
include("utils.jl")
symbolic_has_known_size(x) = !(SU.shape(unwrap(x)) isa SU.Unknown)
abstract type StateMachineOperator end
include("systems/index_cache.jl")
include("systems/abstractsystem.jl")
include("systems/connectors.jl")
include("systems/callbacks.jl")
include("systems/system.jl")
include("systems/codegen_utils.jl")
include("linearization.jl")
include("systems/nonlinear/initializesystem.jl")
include("systems/sparsematrixclil.jl")
include("clock.jl")
include("discretedomain.jl")
include("systems/systemstructure.jl")
include("systems/systems.jl")
include("systems/alias_elimination.jl")
include("structural_transformation/StructuralTransformations.jl")
@reexport using .StructuralTransformations
include("inputoutput.jl")
const t_nounits = let
    only(@independent_variables t)
end
const t = let
    only(@independent_variables t [unit = DQ.u"s"])
end
const D_nounits = Differential(t_nounits)
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
export Sample, Hold, Shift, ShiftIndex, sampletime, SampleTime
export Clock, SolverStepClock, TimeDomain
export MTKParameters, reorder_dimension_by_tunables!, reorder_dimension_by_tunables
export HomotopyContinuationProblem
export AnalysisPoint, get_sensitivity_function, get_comp_sensitivity_function,
       get_looptransfer_function, get_sensitivity, get_comp_sensitivity, get_looptransfer,
       open_loop
function FMIComponent end
const set_scalar_metadata = setmetadata
function __init__()
    SU.hashcons(unwrap(t_nounits), true)
    SU.hashcons(unwrap(t), true)
    SU.hashcons(COMMON_NOTHING, true)
    SU.hashcons(COMMON_MISSING, true)
end
PrecompileTools.@compile_workload begin
     using ModelingToolkit
    @variables x(ModelingToolkit.t_nounits) y(ModelingToolkit.t_nounits)
    isequal(ModelingToolkit.D_nounits.x, ModelingToolkit.t_nounits)
    sys = System([ModelingToolkit.D_nounits(x) ~ x * y, y ~ 3x + 4 * D(y)], ModelingToolkit.t_nounits, [x, y], Num[]; name = :sys)
    TearingState(sys)
end
end
