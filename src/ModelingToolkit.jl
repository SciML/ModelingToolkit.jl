module ModelingToolkit

export Operation, Expression
export calculate_jacobian, generate_jacobian, generate_function
export independent_variables, dependent_variables, parameters
export @register


using DiffEqBase
using StaticArrays, LinearAlgebra

using MacroTools
import MacroTools: splitdef, combinedef

abstract type Expression <: Number end
abstract type AbstractSystem end

Base.promote_rule(::Type{<:Number},::Type{<:Expression}) = Expression
Base.zero(::Type{<:Expression}) = Constant(0)
Base.one(::Type{<:Expression}) = Constant(1)

function calculate_jacobian end
function generate_jacobian end
function generate_function end

function independent_variables end
function dependent_variables   end
function parameters            end

@enum FunctionVersion ArrayFunction=1 SArrayFunction=2

include("variables.jl")
include("operations.jl")
include("differentials.jl")
include("equations.jl")
include("systems/diffeqs/diffeqsystem.jl")
include("systems/diffeqs/first_order_transform.jl")
include("systems/nonlinear/nonlinear_system.jl")
include("function_registration.jl")
include("simplify.jl")
include("utils.jl")

end # module
