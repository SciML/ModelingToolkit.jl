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
                 PeriodicClock, Clock, SolverStepClock, Continuous
using Distributed
import JuliaFormatter
using MLStyle
using NonlinearSolve
import SCCNonlinearSolve
using Reexport
using RecursiveArrayTools
import Graphs: SimpleDiGraph, add_edge!, incidence_matrix
import BlockArrays: BlockArray, BlockedArray, Block, blocksize, blocksizes, blockpush!,
                    undef_blocks, blocks
import CommonSolve
import EnumX

using RuntimeGeneratedFunctions
using RuntimeGeneratedFunctions: drop_expr

using Symbolics: degree
using Symbolics: _parse_vars, value, @derivatives, get_variables,
                 exprs_occur_in, symbolic_linear_solve, build_expr, unwrap, wrap,
                 VariableSource, getname, variable, Connection, connect,
                 NAMESPACE_SEPARATOR, set_scalar_metadata, setdefaultval,
                 initial_state, transition, activeState, entry, hasnode,
                 ticksInState, timeInState, fixpoint_sub, fast_substitute,
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
export independent_variables, unknowns, parameters, full_parameters, continuous_events,
       discrete_events
@reexport using Symbolics
@reexport using UnPack
RuntimeGeneratedFunctions.init(@__MODULE__)

import DynamicQuantities, Unitful
const DQ = DynamicQuantities

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

"""
$(TYPEDEF)

TODO
"""
abstract type AbstractSystem end
abstract type AbstractTimeDependentSystem <: AbstractSystem end
abstract type AbstractTimeIndependentSystem <: AbstractSystem end
abstract type AbstractODESystem <: AbstractTimeDependentSystem end
abstract type AbstractMultivariateSystem <: AbstractSystem end
abstract type AbstractOptimizationSystem <: AbstractTimeIndependentSystem end

function independent_variable end

# this has to be included early to deal with dependency issues
include("structural_transformation/bareiss.jl")
function complete end
function var_derivative! end
function var_derivative_graph! end
include("bipartite_graph.jl")
using .BipartiteGraphs

include("variables.jl")
include("parameters.jl")
include("independent_variables.jl")
include("constants.jl")

include("utils.jl")
include("domains.jl")

include("systems/index_cache.jl")
include("systems/parameter_buffer.jl")
include("systems/abstractsystem.jl")
include("systems/model_parsing.jl")
include("systems/connectors.jl")
include("systems/analysis_points.jl")
include("systems/imperative_affect.jl")
include("systems/callbacks.jl")
include("systems/codegen_utils.jl")
include("systems/problem_utils.jl")

include("systems/nonlinear/nonlinearsystem.jl")
include("systems/nonlinear/homotopy_continuation.jl")
include("systems/diffeqs/odesystem.jl")
include("systems/diffeqs/sdesystem.jl")
include("systems/diffeqs/abstractodesystem.jl")
include("systems/nonlinear/modelingtoolkitize.jl")
include("systems/nonlinear/initializesystem.jl")
include("systems/diffeqs/first_order_transform.jl")
include("systems/diffeqs/modelingtoolkitize.jl")
include("systems/diffeqs/basic_transformations.jl")

include("systems/discrete_system/discrete_system.jl")

include("systems/jumps/jumpsystem.jl")

include("systems/optimization/constraints_system.jl")
include("systems/optimization/optimizationsystem.jl")
include("systems/optimization/modelingtoolkitize.jl")

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

for S in subtypes(ModelingToolkit.AbstractSystem)
    S = nameof(S)
    @eval convert_system(::Type{<:$S}, sys::$S) = sys
end

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

PrecompileTools.@compile_workload begin
    using ModelingToolkit
    @variables x(ModelingToolkit.t_nounits)
    @named sys = ODESystem([ModelingToolkit.D_nounits(x) ~ -x], ModelingToolkit.t_nounits)
    prob = ODEProblem(structural_simplify(sys), [x => 30.0], (0, 100), [], jac = true)
end

export AbstractTimeDependentSystem,
       AbstractTimeIndependentSystem,
       AbstractMultivariateSystem

export ODESystem,
       ODEFunction, ODEFunctionExpr, ODEProblemExpr, convert_system,
       add_accumulations, System
export DAEFunctionExpr, DAEProblemExpr
export SDESystem, SDEFunction, SDEFunctionExpr, SDEProblemExpr
export SystemStructure
export DiscreteSystem, DiscreteProblem, DiscreteFunction, DiscreteFunctionExpr
export JumpSystem
export ODEProblem, SDEProblem
export NonlinearFunction, NonlinearFunctionExpr
export NonlinearProblem, NonlinearProblemExpr
export IntervalNonlinearFunction, IntervalNonlinearFunctionExpr
export IntervalNonlinearProblem, IntervalNonlinearProblemExpr
export OptimizationProblem, OptimizationProblemExpr, constraints
export SteadyStateProblem, SteadyStateProblemExpr
export JumpProblem
export NonlinearSystem, OptimizationSystem, ConstraintsSystem
export alias_elimination, flatten
export connect, domain_connect, @connector, Connection, AnalysisPoint, Flow, Stream,
       instream
export initial_state, transition, activeState, entry, ticksInState, timeInState
export @component, @mtkmodel, @mtkbuild
export isinput, isoutput, getbounds, hasbounds, getguess, hasguess, isdisturbance,
       istunable, getdist, hasdist,
       tunable_parameters, isirreducible, getdescription, hasdescription,
       hasunit, getunit, hasconnect, getconnect,
       hasmisc, getmisc, state_priority
export ode_order_lowering, dae_order_lowering, liouville_transform
export PDESystem
export Differential, expand_derivatives, @derivatives
export Equation, ConstrainedEquation
export Term, Sym
export SymScope, LocalScope, ParentScope, DelayParentScope, GlobalScope
export independent_variable, equations, controls, observed, full_equations
export initialization_equations, guesses, defaults, parameter_dependencies, hierarchy
export structural_simplify, expand_connections, linearize, linearization_function

export calculate_jacobian, generate_jacobian, generate_function, generate_custom_function
export calculate_control_jacobian, generate_control_jacobian
export calculate_tgrad, generate_tgrad
export calculate_gradient, generate_gradient
export calculate_factorized_W, generate_factorized_W
export calculate_hessian, generate_hessian
export calculate_massmatrix, generate_diffusion_function
export stochastic_integral_transform
export TearingState

export BipartiteGraph, equation_dependencies, variable_dependencies
export eqeq_dependencies, varvar_dependencies
export asgraph, asdigraph

export toexpr, get_variables
export simplify, substitute
export build_function
export modelingtoolkitize
export generate_initializesystem

export alg_equations, diff_equations, has_alg_equations, has_diff_equations
export get_alg_eqs, get_diff_eqs, has_alg_eqs, has_diff_eqs

export @variables, @parameters, @independent_variables, @constants, @brownian
export @named, @nonamespace, @namespace, extend, compose, complete
export debug_system

#export Continuous, Discrete, sampletime, input_timedomain, output_timedomain
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

end # module
