"""
    $(TYPEDEF)

Wrapper over an `AbstractDict{K, SymbolicT} where {K}` which matches higher order
derivatives if the first-order derivative is present in the wrapped dictionary.
Specifically, if `Differential(t, 1)(x)` is present in the wrapped dictionary with value
`v`, then `Differential(t, n::Int)(x)` maps to `Differential(t, n - 1)(v)`. This only
affects `get`, `getindex` and `haskey`. All other methods fall back to the wrapped
dictionary.
"""
struct DerivativeDict{K, D <: AbstractDict{K, SymbolicT}} <: AbstractDict{K, SymbolicT}
    dict::D
end

DerivativeDict() = DerivativeDict(Dict{SymbolicT, SymbolicT}())
DerivativeDict(args::Pair...) = DerivativeDict(Dict(args...))
DerivativeDict{K}(args...) where {K} = DerivativeDict(Dict{K, SymbolicT}(args...))

Base.copy(dd::DerivativeDict) = DerivativeDict(copy(dd.dict))
function Base.empty(dd::DerivativeDict, ::Type{K}, ::Type{V}) where {K, V}
    DerivativeDict(empty(dd.dict, K, V))
end

struct __DDSentinel end
const DD_SENTINEL = __DDSentinel()

function Base.get(def::Base.Callable, dd::DerivativeDict, k::SymbolicT)
    res = get(dd.dict, k, DD_SENTINEL)
    res === DD_SENTINEL || return res
    Moshi.Match.@match k begin
        BSImpl.Term(; f, args) && if f isa Differential && f.order::Int > 1 end => begin
            order = f.order::Int
            res = get(dd.dict, Differential(f.x, 1)(args[1]), DD_SENTINEL)
            res === DD_SENTINEL && return def()
            res = res::SymbolicT
            return Differential(f.x, order - 1)(res)
        end
        _ => return def()
    end
end
function Base.get(f::Base.Callable, dd::DerivativeDict, k::Union{Num, Arr, CallAndWrap})
    return get(f, dd, unwrap(k))
end
Base.get(f::Base.Callable, dd::DerivativeDict, k) = get(f, dd.dict, k)

Base.get(dd::DerivativeDict, k, default) = get(Returns(default), dd, k)

Base.haskey(dd::DerivativeDict, k) = get(Returns(DD_SENTINEL), dd, k) !== DD_SENTINEL

function Base.getindex(dd::DerivativeDict, k)
    res = get(Returns(DD_SENTINEL), dd, k)
    res === DD_SENTINEL && throw(KeyError(k))
    return res::SymbolicT
end

Base.setindex!(dd::DerivativeDict, v, k) = setindex!(dd.dict, unwrap(v), unwrap(k))

Base.isempty(dd::DerivativeDict) = isempty(dd.dict)
Base.length(dd::DerivativeDict) = length(dd.dict)
Base.iterate(dd::DerivativeDict, args...) = Base.iterate(dd.dict, args...)
Base.sizehint!(dd::DerivativeDict, n; kw...) = sizehint!(dd.dict, n; kw...)
Base.empty!(dd::DerivativeDict) = empty!(dd.dict)
