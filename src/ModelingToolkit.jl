module ModelingToolkit

using DiffEqBase
using StaticArrays, LinearAlgebra

using MacroTools
import MacroTools: splitdef, combinedef

abstract type Expression <: Number end
abstract type AbstractOperation <: Expression end
abstract type AbstractComponent <: Expression end

include("variables.jl")

Base.promote_rule(::Type{T},::Type{T2}) where {T<:Number,T2<:Expression} = Expression
Base.zero(::Type{<:Expression}) = Constant(0)
Base.one(::Type{<:Expression}) = Constant(1)
Base.convert(::Type{Variable},x::Int64) = Constant(x)

function caclulate_jacobian end

@enum FunctionVersion ArrayFunction=1 SArrayFunction=2

include("operations.jl")
include("differentials.jl")
include("equations.jl")
include("systems/systems.jl")
include("systems/diffeqs/diffeqsystem.jl")
include("systems/diffeqs/first_order_transform.jl")
include("systems/nonlinear/nonlinear_system.jl")
include("function_registration.jl")
include("simplify.jl")
include("utils.jl")

export Operation, Expression, AbstractComponent, AbstractDomain
export @register
end # module
