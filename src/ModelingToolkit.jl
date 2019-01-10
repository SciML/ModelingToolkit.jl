module ModelingToolkit

using DiffEqBase
using StaticArrays, LinearAlgebra
import Terms: Term, root, children, is_branch, @term

using MacroTools

abstract type Expression <: Number end
abstract type AbstractComponent <: Expression end
abstract type AbstractSystem end
abstract type AbstractDomain end

include("domains.jl")
include("variables.jl")

Base.promote_rule(::Type{T},::Type{T2}) where {T<:Number,T2<:Expression} = Expression
Base.zero(::Type{<:Expression}) = @term(0)
Base.one(::Type{<:Expression}) = @term(1)
Base.convert(::Type{Variable}, x::Int64) = convert(Term, x)

function caclulate_jacobian end

@enum FunctionVersions ArrayFunction=1 SArrayFunction=2

include("differentials.jl")
include("systems/diffeqs/diffeqsystem.jl")
include("systems/diffeqs/first_order_transform.jl")
include("systems/nonlinear/nonlinear_system.jl")
include("simplify.jl")
include("utils.jl")

export Term, @term, Expression, AbstractComponent, AbstractDomain
end # module
