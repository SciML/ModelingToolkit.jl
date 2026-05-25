module MTKBaseEnzymeExt

using ModelingToolkitBase: GeneratedFunctionWrapper, ObservedFunctionCache,
    MissingGuessValue, System
import Enzyme: EnzymeRules

# Mirror of the declarations in `MTKMooncakeExt` that mark non-differentiable
# structural / configuration types as inactive. Mooncake handles this through
# `Mooncake.tangent_type(::Type{T}) = NoTangent`; Enzyme uses
# `EnzymeRules.inactive_type`.
#
# Without these, Enzyme's runtime-activity dispatch and `make_zero` traversal
# recurse into self-referential or otherwise-unhandleable fields when a
# closure (or problem/solution struct) transitively captures one of these
# values, producing activity-analysis errors or `MethodError`s in
# `create_activity_wrapper`.

# `AbstractSystem` is already declared `inactive_type` in
# `ModelingToolkitBase.jl` itself (EnzymeCore is a hard dep there); this
# subtype declaration is redundant but explicit, and keeps the extension in
# direct one-to-one parity with `MTKMooncakeExt`.
EnzymeRules.inactive_type(::Type{<:System}) = true

# `GeneratedFunctionWrapper` wraps RuntimeGeneratedFunctions - the model RHS
# closure, not differentiable data.
EnzymeRules.inactive_type(::Type{<:GeneratedFunctionWrapper}) = true

# `ObservedFunctionCache` contains a `System` reference and a `Dict{Any,Any}`
# of cached observed-function closures. Not a derivative carrier.
EnzymeRules.inactive_type(::Type{<:ObservedFunctionCache}) = true

# `MissingGuessValue.Type` is a Moshi `@data` tagged union
# (`Constant{<:Number} | Random{<:AbstractRNG} | Error`) used as a
# configuration enum for initialization. Not differentiable.
EnzymeRules.inactive_type(::Type{<:MissingGuessValue.Type}) = true

end
