module MTKFunctorsExt

using Functors: Functors
using ModelingToolkitBase: NonNumericWrapper

# `NonNumericWrapper` exists so that AD backends can treat the nonnumeric portion of
# `MTKParameters` as inert. SciMLSensitivity walks a parameter object together with its
# cotangent through `Functors.fmap`, and the cotangent of this slot is `nothing`. Functors
# destructures every argument with the functor of the *first* one, so it needs to know how
# to take `(buffer = ...)` children off `nothing` and off cotangent carriers such as
# `NamedTuple`s or `Tangent`s, not just off the wrapper itself. The single-argument walks
# (`allocate_vjp`, `allocate_zeros`, ...) rely on the wrapper having children, so it must
# not be declared a leaf: they would otherwise recurse forever on it.
Functors.functor(::Type{<:NonNumericWrapper}, x::NonNumericWrapper) =
    (; buffer = x.buffer), y -> NonNumericWrapper(y.buffer)
Functors.functor(::Type{<:NonNumericWrapper}, ::Nothing) = (; buffer = nothing), identity
Functors.functor(::Type{<:NonNumericWrapper}, x) = (; buffer = x.buffer), identity

end
