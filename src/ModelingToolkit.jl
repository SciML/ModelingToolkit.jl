module ModelingToolkit

using DiffEqBase, Distributed
using StaticArrays, LinearAlgebra, SparseArrays
using Latexify
using MacroTools

using MacroTools
import MacroTools: splitdef, combinedef, postwalk, striplines
import GeneralizedGenerated
using DocStringExtensions
using Base: RefValue

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

Calculate the jacobian matrix of a system.

Returns a matrix of [`Expression`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_jacobian end

"""
$(TYPEDSIGNATURES)

Generate a function to calculate the Jacobian of the system.
"""
function generate_jacobian end

"""
$(TYPEDSIGNATURES)

Generate a function to evaluate the system's equations.
"""
function generate_function end

"""
$(TYPEDSIGNATURES)

Get the set of independent variables for the given system.
"""
function independent_variables end

"""
$(TYPEDSIGNATURES)

Get the set of dependent variables for the given system.
"""
function dependent_variables end

"""
$(TYPEDSIGNATURES)

Get the set of parameters variables for the given system.
"""
function parameters end

include("variables.jl")
include("operations.jl")
include("differentials.jl")
include("equations.jl")
include("function_registration.jl")
include("simplify.jl")
include("utils.jl")
include("direct.jl")
include("domains.jl")
include("systems/diffeqs/odesystem.jl")
include("systems/diffeqs/sdesystem.jl")
include("systems/diffeqs/abstractodesystem.jl")
include("systems/diffeqs/first_order_transform.jl")
include("systems/diffeqs/modelingtoolkitize.jl")
include("systems/nonlinear/nonlinear_system.jl")
include("systems/pde/pdesystem.jl")
include("systems/reaction/reactionsystem.jl")
include("latexify_recipes.jl")
include("build_function.jl")

export ODESystem, ODEFunction
export SDESystem, SDEFunction
export NonlinearSystem
export ode_order_lowering
export PDESystem
export Reaction, ReactionSystem
export Differential, expand_derivatives, @derivatives
export IntervalDomain, ProductDomain, ⊗, CircleDomain
export Equation, ConstrainedEquation
export simplify_constants

export Operation, Expression
export calculate_jacobian, generate_jacobian, generate_function
export calculate_massmatrix, generate_diffusion_function
export independent_variable, states, parameters, equations
export simplified_expr, eval_function
export @register, @I
export modelingtoolkitize
export Variable, @variables, @parameters
end # module
