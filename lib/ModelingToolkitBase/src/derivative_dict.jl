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

"""
    $(TYPEDEF)

Identical to `DerivativeDict`, but `Differential(t, n)(x)` maps to
`expand_derivatives(Differential(t, n-1)(v))`. Mutates the wrapped dictionary in `get`
for efficiency - the `expand_derivatives` result is cached. This also means that `getindex`
mutates.
"""
struct ExpandDerivativeDict{K, D <: AbstractDict{K, SymbolicT}} <: AbstractDict{K, SymbolicT}
    dict::D
end

struct __DDSentinel end
const DD_SENTINEL = __DDSentinel()

for T in [:DerivativeDict, :ExpandDerivativeDict]
    @eval begin
        $T() = $T(Dict{SymbolicT, SymbolicT}())
        $T(args::Pair...) = $T(Dict(args...))
        $T{K}(args...) where {K} = $T(Dict{K, SymbolicT}(args...))

        Base.copy(dd::$T) = $T(copy(dd.dict))
        function Base.empty(dd::$T, ::Type{K}, ::Type{V}) where {K, V}
            return $T(empty(dd.dict, K, V))
        end

        Base.haskey(dd::$T, k) = get(Returns(DD_SENTINEL), dd, k) !== DD_SENTINEL

        Base.setindex!(dd::$T, v, k) = setindex!(dd.dict, unwrap(v), unwrap(k))

        Base.isempty(dd::$T) = isempty(dd.dict)
        Base.length(dd::$T) = length(dd.dict)
        Base.iterate(dd::$T, args...) = Base.iterate(dd.dict, args...)
        Base.sizehint!(dd::$T, n; kw...) = sizehint!(dd.dict, n; kw...)
        Base.empty!(dd::$T) = empty!(dd.dict)

        function Base.get(def::Base.Callable, dd::$T, k::SymbolicT)
            res = get(dd.dict, k, DD_SENTINEL)
            res === DD_SENTINEL || return res
            return Moshi.Match.@match k begin
                BSImpl.Term(; f, args) && if f isa Differential && f.order::Int > 1 end => begin
                    order = f.order::Int
                    res = get(dd.dict, Differential(f.x, 1)(args[1]), DD_SENTINEL)
                    res === DD_SENTINEL && return def()
                    res = res::SymbolicT
                    $(
                        if T === :DerivativeDict
                            :(return Differential(f.x, order - 1)(res))
                        elseif T === :ExpandDerivativeDict
                            quote
                                res = expand_derivatives(Differential(f.x, order - 1)(res))
                                dd.dict[k] = res
                                return res
                            end
                        else
                            error("Not handled")
                        end
                    )
                end
                _ => return def()
            end
        end
        function Base.get(f::Base.Callable, dd::$T, k::Union{Num, Arr, CallAndWrap})
            return get(f, dd, unwrap(k))
        end
        Base.get(f::Base.Callable, dd::$T, k) = get(f, dd.dict, k)

        Base.get(dd::$T, k, default) = get(Returns(default), dd, k)

        function Base.getindex(dd::$T, k)
            res = get(Returns(DD_SENTINEL), dd, k)
            res === DD_SENTINEL && throw(KeyError(k))
            return res::SymbolicT
        end

    end
end

