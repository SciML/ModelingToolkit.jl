"""
$(DocStringExtensions.README)
"""
module ModelingToolkit
using DocStringExtensions
using AbstractTrees
using DiffEqBase, SciMLBase, ForwardDiff, Reexport
using Distributed
using StaticArrays, LinearAlgebra, SparseArrays, LabelledArrays
using InteractiveUtils
using Latexify, Unitful, ArrayInterfaceCore
using MacroTools
@reexport using UnPack
using Setfield, ConstructionBase
using JumpProcesses
using DataStructures
using SpecialFunctions, NaNMath
using RuntimeGeneratedFunctions
using Base.Threads
using DiffEqCallbacks
using Graphs
import MacroTools: splitdef, combinedef, postwalk, striplines
import Libdl
using DocStringExtensions
using Base: RefValue
using Combinatorics
import IfElse
import Distributions
import FunctionWrappersWrappers

RuntimeGeneratedFunctions.init(@__MODULE__)

using RecursiveArrayTools

import SymbolicUtils
import SymbolicUtils: istree, arguments, operation, similarterm, promote_symtype,
                      Symbolic, Term, Add, Mul, Pow, Sym, FnType,
                      @rule, Rewriters, substitute, metadata
using SymbolicUtils.Code
import SymbolicUtils.Code: toexpr
import SymbolicUtils.Rewriters: Chain, Postwalk, Prewalk, Fixpoint
import JuliaFormatter

using Reexport
using Symbolics: degree
@reexport using Symbolics
export @derivatives
using Symbolics: _parse_vars, value, @derivatives, get_variables,
                 exprs_occur_in, solve_for, build_expr, unwrap, wrap,
                 VariableSource, getname, variable, Connection, connect,
                 NAMESPACE_SEPARATOR
import Symbolics: rename, get_variables!, _solve, hessian_sparsity,
                  jacobian_sparsity, isaffine, islinear, _iszero, _isone,
                  tosymbol, lower_varname, diff2term, var_from_nested_derivative,
                  BuildTargets, JuliaTarget, StanTarget, CTarget, MATLABTarget,
                  ParallelForm, SerialForm, MultithreadedForm, build_function,
                  rhss, lhss, prettify_expr, gradient,
                  jacobian, hessian, derivative, sparsejacobian, sparsehessian,
                  substituter, scalarize, getparent

import DiffEqBase: @add_kwonly

import Graphs: SimpleDiGraph, add_edge!, incidence_matrix

for fun in [:toexpr]
    @eval begin
        function $fun(eq::Equation; kw...)
            Expr(:(=), $fun(eq.lhs; kw...), $fun(eq.rhs; kw...))
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

"""
$(TYPEDSIGNATURES)

Get the set of independent variables for the given system.
"""
function independent_variables end

function independent_variable end

"""
$(TYPEDSIGNATURES)

Get the set of states for the given system.
"""
function states end

"""
$(TYPEDSIGNATURES)

Get the set of parameters variables for the given system.
"""
function parameters end

# this has to be included early to deal with depency issues
include("structural_transformation/bareiss.jl")
include("bipartite_graph.jl")
using .BipartiteGraphs

include("variables.jl")
include("parameters.jl")

include("utils.jl")
include("domains.jl")

include("systems/abstractsystem.jl")
include("systems/connectors.jl")
include("systems/callbacks.jl")

include("systems/diffeqs/odesystem.jl")
include("systems/diffeqs/sdesystem.jl")
include("systems/diffeqs/abstractodesystem.jl")
include("systems/diffeqs/first_order_transform.jl")
include("systems/diffeqs/modelingtoolkitize.jl")
include("systems/diffeqs/basic_transformations.jl")

include("systems/jumps/jumpsystem.jl")

include("systems/nonlinear/nonlinearsystem.jl")
include("systems/nonlinear/modelingtoolkitize.jl")

include("systems/optimization/optimizationsystem.jl")

include("systems/pde/pdesystem.jl")

include("systems/sparsematrixclil.jl")
include("systems/discrete_system/discrete_system.jl")
include("systems/validation.jl")
include("systems/dependency_graphs.jl")
include("systems/systemstructure.jl")
using .SystemStructures

include("systems/alias_elimination.jl")
include("structural_transformation/StructuralTransformations.jl")

@reexport using .StructuralTransformations
include("inputoutput.jl")

for S in subtypes(ModelingToolkit.AbstractSystem)
    S = nameof(S)
    @eval convert_system(::Type{<:$S}, sys::$S) = sys
end

export AbstractTimeDependentSystem, AbstractTimeIndependentSystem,
       AbstractMultivariateSystem
export ODESystem, ODEFunction, ODEFunctionExpr, ODEProblemExpr, convert_system
export DAEFunctionExpr, DAEProblemExpr
export SDESystem, SDEFunction, SDEFunctionExpr, SDEProblemExpr
export SystemStructure
export JumpSystem
export ODEProblem, SDEProblem
export NonlinearFunction, NonlinearFunctionExpr
export NonlinearProblem, BlockNonlinearProblem, NonlinearProblemExpr
export OptimizationProblem, OptimizationProblemExpr, constraints
export AutoModelingToolkit
export SteadyStateProblem, SteadyStateProblemExpr
export JumpProblem, DiscreteProblem
export NonlinearSystem, OptimizationSystem
export alias_elimination, flatten
export connect, @connector, Connection, Flow, Stream, instream
export isinput, isoutput, getbounds, hasbounds, isdisturbance, istunable, getdist, hasdist,
       tunable_parameters, isirreducible, getdescription, hasdescription, isbinaryvar,
       isintegervar
export ode_order_lowering, dae_order_lowering, liouville_transform
export PDESystem
export Differential, expand_derivatives, @derivatives
export Equation, ConstrainedEquation
export Term, Sym
export SymScope, LocalScope, ParentScope, GlobalScope
export independent_variables, independent_variable, states, parameters, equations, controls,
       observed, structure, full_equations
export structural_simplify, expand_connections, linearize, linear_statespace
export DiscreteSystem, DiscreteProblem

export calculate_jacobian, generate_jacobian, generate_function
export calculate_control_jacobian, generate_control_jacobian
export calculate_tgrad, generate_tgrad
export calculate_gradient, generate_gradient
export calculate_factorized_W, generate_factorized_W
export calculate_hessian, generate_hessian
export calculate_massmatrix, generate_diffusion_function
export stochastic_integral_transform
export TearingState, StateSelectionState
export generate_difference_cb

export BipartiteGraph, equation_dependencies, variable_dependencies
export eqeq_dependencies, varvar_dependencies
export asgraph, asdigraph

export toexpr, get_variables
export simplify, substitute
export build_function
export modelingtoolkitize
export @variables, @parameters
export @named, @nonamespace, @namespace, extend, compose

end # module
