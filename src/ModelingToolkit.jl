module ModelingToolkit

export Operation, Expression
export calculate_jacobian, generate_jacobian, generate_function
export independent_variables, dependent_variables, parameters
export simplified_expr, eval_function
export @register, @I
export modelingtoolkitize


using DiffEqBase, Distributed
using StaticArrays, LinearAlgebra

using MacroTools
import MacroTools: splitdef, combinedef
import GeneralizedGenerated
using DocStringExtensions

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

@enum FunctionVersion ArrayFunction=1 SArrayFunction=2

include("variables.jl")
include("operations.jl")
include("differentials.jl")
include("equations.jl")
include("function_registration.jl")
include("simplify.jl")
include("utils.jl")
include("direct.jl")
include("domains.jl")
include("systems/diffeqs/diffeqsystem.jl")
include("systems/diffeqs/first_order_transform.jl")
include("systems/nonlinear/nonlinear_system.jl")
include("systems/pde/pdesystem.jl")

end # module
