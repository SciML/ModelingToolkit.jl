"""
$(DocStringExtensions.README)
"""
module ModelingToolkit
using PrecompileTools, Reexport
@recompile_invalidations begin
    using DocStringExtensions
    using Compat
    using AbstractTrees
    using DiffEqBase, SciMLBase, ForwardDiff
    using SciMLBase: StandardODEProblem, StandardNonlinearProblem, handle_varmap
    using Distributed
    using StaticArrays, LinearAlgebra, SparseArrays, LabelledArrays
    using InteractiveUtils
    using Latexify, Unitful, ArrayInterface
    using MacroTools
    using Setfield, ConstructionBase
    using JumpProcesses
    using DataStructures
    using SpecialFunctions, NaNMath
    using RuntimeGeneratedFunctions
    using RuntimeGeneratedFunctions: drop_expr
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
    using URIs: URI

    using RecursiveArrayTools

    using SymbolicIndexingInterface
    export independent_variables, states, parameters
    import SymbolicUtils
    import SymbolicUtils: istree, arguments, operation, similarterm, promote_symtype,
        Symbolic, isadd, ismul, ispow, issym, FnType,
        @rule, Rewriters, substitute, metadata, BasicSymbolic,
        Sym, Term
    using SymbolicUtils.Code
    import SymbolicUtils.Code: toexpr
    import SymbolicUtils.Rewriters: Chain, Postwalk, Prewalk, Fixpoint
    import JuliaFormatter

    using MLStyle

    using Reexport
    using Symbolics
    using Symbolics: degree
    using Symbolics: _parse_vars, value, @derivatives, get_variables,
        exprs_occur_in, solve_for, build_expr, unwrap, wrap,
        VariableSource, getname, variable, Connection, connect,
        NAMESPACE_SEPARATOR, set_scalar_metadata, setdefaultval
    import Symbolics: rename, get_variables!, _solve, hessian_sparsity,
        jacobian_sparsity, isaffine, islinear, _iszero, _isone,
        tosymbol, lower_varname, diff2term, var_from_nested_derivative,
        BuildTargets, JuliaTarget, StanTarget, CTarget, MATLABTarget,
        ParallelForm, SerialForm, MultithreadedForm, build_function,
        rhss, lhss, prettify_expr, gradient,
        jacobian, hessian, derivative, sparsejacobian, sparsehessian,
        substituter, scalarize, getparent

    import DiffEqBase: @add_kwonly
    import OrdinaryDiffEq

    import Graphs: SimpleDiGraph, add_edge!, incidence_matrix
end

@reexport using Symbolics
@reexport using UnPack
RuntimeGeneratedFunctions.init(@__MODULE__)

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
include("constants.jl")

include("utils.jl")
include("domains.jl")

include("systems/abstractsystem.jl")
include("systems/model_parsing.jl")
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
include("systems/nonlinear/initializesystem.jl")

include("systems/optimization/constraints_system.jl")
include("systems/optimization/optimizationsystem.jl")
include("systems/optimization/modelingtoolkitize.jl")

include("systems/pde/pdesystem.jl")

include("systems/sparsematrixclil.jl")
include("systems/discrete_system/discrete_system.jl")
include("systems/validation.jl")
include("systems/dependency_graphs.jl")
include("clock.jl")
include("discretedomain.jl")
include("systems/systemstructure.jl")
include("systems/clock_inference.jl")
include("systems/systems.jl")

include("debugging.jl")
include("systems/alias_elimination.jl")
include("structural_transformation/StructuralTransformations.jl")

@reexport using .StructuralTransformations
include("inputoutput.jl")

for S in subtypes(ModelingToolkit.AbstractSystem)
    S = nameof(S)
    @eval convert_system(::Type{<:$S}, sys::$S) = sys
end

PrecompileTools.@compile_workload begin
    using ModelingToolkit
    @variables t x(t)
    D = Differential(t)
    @named sys = ODESystem([D(x) ~ -x])
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
export JumpSystem
export ODEProblem, SDEProblem
export NonlinearFunction, NonlinearFunctionExpr
export NonlinearProblem, BlockNonlinearProblem, NonlinearProblemExpr
export OptimizationProblem, OptimizationProblemExpr, constraints
export AutoModelingToolkit
export SteadyStateProblem, SteadyStateProblemExpr
export JumpProblem, DiscreteProblem
export NonlinearSystem, OptimizationSystem, ConstraintsSystem
export alias_elimination, flatten
export connect, domain_connect, @connector, Connection, Flow, Stream, instream
export @component, @mtkmodel, @mtkbuild
export isinput, isoutput, getbounds, hasbounds, getguess, hasguess, isdisturbance,
    istunable, getdist, hasdist,
    tunable_parameters, isirreducible, getdescription, hasdescription, isbinaryvar,
    isintegervar
export ode_order_lowering, dae_order_lowering, liouville_transform
export PDESystem
export Differential, expand_derivatives, @derivatives
export Equation, ConstrainedEquation
export Term, Sym
export SymScope, LocalScope, ParentScope, DelayParentScope, GlobalScope
export independent_variable, equations, controls,
    observed, structure, full_equations
export structural_simplify, expand_connections, linearize, linearization_function
export DiscreteSystem,
    DiscreteProblem, DiscreteProblemExpr, DiscreteFunction,
    DiscreteFunctionExpr

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
export initializesystem

export @variables, @parameters, @constants, @brownian
export @named, @nonamespace, @namespace, extend, compose, complete
export debug_system

export show_with_compare

#export Continuous, Discrete, sampletime, input_timedomain, output_timedomain
#export has_discrete_domain, has_continuous_domain
#export is_discrete_domain, is_continuous_domain, is_hybrid_domain
export Sample, Hold, Shift, ShiftIndex
export Clock #, InferredDiscrete,

end # module
