module ModelingToolkit

using DiffEqBase, Distributed
using StaticArrays, LinearAlgebra, SparseArrays
using Latexify, Unitful, ArrayInterface
using MacroTools
using UnPack: @unpack
using DiffEqJump

using Base.Threads
import MacroTools: splitdef, combinedef, postwalk, striplines
import GeneralizedGenerated
using DocStringExtensions
using Base: RefValue

using RecursiveArrayTools

import SymbolicUtils
import SymbolicUtils: to_symbolic, FnType

import TreeViews

"""
$(TYPEDEF)

Base type for a symbolic expression.
"""
abstract type Expression <: Number end

"""
$(TYPEDEF)

TODO
"""
abstract type AbstractSystem end
abstract type AbstractODESystem <: AbstractSystem end

Base.promote_rule(::Type{<:Number},::Type{<:Expression}) = Expression
Base.zero(::Type{<:Expression}) = Constant(0)
Base.one(::Type{<:Expression}) = Constant(1)

"""
$(TYPEDSIGNATURES)

Get the set of independent variables for the given system.
"""
function independent_variables end

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

include("variables.jl")
include("context_dsl.jl")
include("operations.jl")
include("differentials.jl")

function Base.convert(::Type{Variable},x::Operation)
    if x.op isa Variable
        x.op
    elseif x.op isa Differential
        var = x.args[1].op
        rename(var,Symbol(var.name,:ˍ,x.args[1].args[1].op.name))
    else
        throw(error("This Operation is not a Variable"))
    end
end

include("equations.jl")
include("function_registration.jl")
include("simplify.jl")
include("utils.jl")
include("direct.jl")
include("domains.jl")

include("systems/abstractsystem.jl")

include("systems/diffeqs/odesystem.jl")
include("systems/diffeqs/sdesystem.jl")
include("systems/diffeqs/abstractodesystem.jl")
include("systems/diffeqs/first_order_transform.jl")
include("systems/diffeqs/modelingtoolkitize.jl")
include("systems/diffeqs/validation.jl")

include("systems/jumps/jumpsystem.jl")

include("systems/nonlinear/nonlinearsystem.jl")

include("systems/optimization/optimizationsystem.jl")

include("systems/pde/pdesystem.jl")

include("systems/reaction/reactionsystem.jl")
include("systems/dependency_graphs.jl")

include("latexify_recipes.jl")
include("build_function.jl")

export ODESystem, ODEFunction
export SDESystem, SDEFunction
export JumpSystem
export ODEProblem, SDEProblem, NonlinearProblem, OptimizationProblem
export JumpProblem, DiscreteProblem
export NonlinearSystem, OptimizationSystem
export ode_order_lowering
export PDESystem
export Reaction, ReactionSystem
export Differential, expand_derivatives, @derivatives
export IntervalDomain, ProductDomain, ⊗, CircleDomain
export Equation, ConstrainedEquation
export Operation, Expression, Variable
export independent_variable, states, parameters, equations 

export calculate_jacobian, generate_jacobian, generate_function
export calculate_tgrad, generate_tgrad
export calculate_gradient, generate_gradient
export calculate_factorized_W, generate_factorized_W
export calculate_hessian, generate_hessian
export calculate_massmatrix, generate_diffusion_function

export simplified_expr, rename, get_variables
export simplify, substitute
export build_function
export @register
export modelingtoolkitize
export @variables, @parameters
end # module
