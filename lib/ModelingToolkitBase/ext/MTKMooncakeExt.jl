module MTKMooncakeExt

using ModelingToolkitBase: MTKParameters, ParameterIndex, NONNUMERIC_PORTION, AbstractSystem,
    SetInitialUnknowns, GeneratedFunctionWrapper, ObservedFunctionCache, MissingGuessValue
using ModelingToolkitBase: System
import ModelingToolkitBase as MTK
import Mooncake
using Mooncake: @from_rrule, @zero_adjoint, MinimalCtx, NoTangent, NoRData, NoFData, increment_and_get_rdata!
import SymbolicIndexingInterface: remake_buffer

# Mooncake @from_rrule interface for MTKParameters.
function Mooncake.increment_and_get_rdata!(f::Mooncake.FData{@NamedTuple{tunable::Vector{T}, initials::Vector{P}, discrete::NoFData, constant::NoFData, nonnumeric::NoFData, caches::NoFData}}, r::Mooncake.NoRData, t::@NamedTuple{tunable::Vector{T}, initials::Nothing, discrete::Nothing, constant::Nothing, nonnumeric::Nothing, caches::Nothing}) where {T <: Base.IEEEFloat, P <: Base.IEEEFloat}
    Mooncake.increment_and_get_rdata!(f.data.tunable, r, t.tunable)
    Mooncake.increment_and_get_rdata!(f.data.caches, r, t.caches)
    return NoRData()
end

# Fix StackOverflow in tangent_type for recursive/structural types.
# These types are model structure descriptors, not differentiable quantities.

# ImmutableDict has a self-referential `parent` field that causes infinite recursion.
Mooncake.tangent_type(::Type{<:Base.ImmutableDict}) = NoTangent

# System has self-referential fields: systems::Vector{System},
# parent::Union{Nothing, System}, initializesystem::Union{Nothing, System}
Mooncake.tangent_type(::Type{<:System}) = NoTangent

# GeneratedFunctionWrapper wraps RuntimeGeneratedFunctions - model RHS, not differentiable.
Mooncake.tangent_type(::Type{<:GeneratedFunctionWrapper}) = NoTangent

# ObservedFunctionCache contains System reference and Dict{Any,Any}.
Mooncake.tangent_type(::Type{<:ObservedFunctionCache}) = NoTangent

# MissingGuessValue is a Moshi @data tagged union with a
# Union{Error, Constant, Random} storage field. Mooncake's FData/RData
# decomposition cannot handle Union-typed tangent fields where the member
# tangent types are distinct Tangent{NamedTuple{...}} structs.
# This is not differentiable anyway — it's a configuration enum.
Mooncake.tangent_type(::Type{<:MissingGuessValue.Type}) = NoTangent

# Port ChainRules rrules to Mooncake using @from_rrule.
# These mirror the rules in MTKChainRulesCoreExt.jl.

# MTKParameters constructor: rrule(::Type{MTKParameters}, tunables, args...)
# The constructor is called with 6 args (tunables, initials, discrete, constant, nonnumeric, caches)
@from_rrule(
    MinimalCtx,
    Tuple{Type{MTKParameters}, Any, Any, Any, Any, Any, Any},
    true,
)

# remake_buffer with MTKParameters
@from_rrule(
    MinimalCtx,
    Tuple{typeof(remake_buffer), Any, MTKParameters, Any, Any},
    true,
)

# SetInitialUnknowns callable struct
@from_rrule(
    MinimalCtx,
    Tuple{SetInitialUnknowns, MTKParameters, Any},
    true,
)
@from_rrule(
    MinimalCtx,
    Tuple{SetInitialUnknowns, AbstractVector, Any},
    true,
)

# Base.getproperty on AbstractSystem is not differentiable
@zero_adjoint MinimalCtx Tuple{typeof(Base.getproperty), AbstractSystem, Symbol}

end
