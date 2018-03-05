module SciCompDSL

using DiffEqBase
import MacroTools: splitdef, combinedef
import IterTools: product

abstract type Expression <: Number end
abstract type AbstractOperation <: Expression end
abstract type AbstractOperator <: Expression end
abstract type AbstractSystem end

include("variables.jl")

Base.promote_rule(::Type{T},::Type{T2}) where {T<:Number,T2<:Expression} = Expression
Base.one(::Type{T}) where T<:Expression = Constant(1)
Base.zero(::Type{T}) where T<:Expression = Constant(0)

include("operations.jl")
include("operators.jl")
include("systems.jl")
include("function_registration.jl")
include("simplify.jl")
include("utils.jl")

export Operation, Expression, AbstractOperator
export @register
end # module
