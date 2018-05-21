module ModelingToolkit

using DiffEqBase
using StaticArrays
using Requires

import MacroTools: splitdef, combinedef
import IterTools: product

@require Reduce begin
    using Reduce
    using ReduceLinAlg
    import Reduce: RExpr
    Reduce.Rational(false)
    Algebra.operator(:\,:inv,:identity)
    Algebra.rlet(:(identity(~x))=>:x)
    Algebra.rlet(:(inv(~x))=>:(x^(-1)))
end

abstract type Expression <: Number end
abstract type AbstractOperation <: Expression end
abstract type AbstractOperator <: Expression end
abstract type AbstractComponent <: Expression end
abstract type AbstractConnection <: Expression end
abstract type AbstractSystem end
abstract type AbstractDomain end

include("domains.jl")
include("variables.jl")

Base.promote_rule(::Type{T},::Type{T2}) where {T<:Number,T2<:Expression} = Expression
Base.one(::Type{T}) where T<:Expression = Constant(1)
Base.zero(::Type{T}) where T<:Expression = Constant(0)
Base.convert(::Type{Variable},x::Int64) = Constant(x)

function caclulate_jacobian end

include("operations.jl")
include("operators.jl")
include("systems/diffeqs/diffeqsystem.jl")
include("systems/diffeqs/first_order_transform.jl")
include("systems/nonlinear/nonlinear_system.jl")
include("function_registration.jl")
include("connections.jl")
include("simplify.jl")
include("utils.jl")

export Operation, Expression, AbstractOperator, AbstractComponent, AbstractDomain
export @register
end # module
