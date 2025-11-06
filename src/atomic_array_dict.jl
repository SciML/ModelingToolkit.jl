"""
    $(TYPEDEF)

Wrapper over an `AbstractDict{SymbolicT, V} where {V}` which disallows keys that are
indexed array symbolics. Specifically, if `@variables x[1:4]` exists, then `x` can be
a key but `x[1]` cannot.
"""
struct AtomicArrayDict{V, D <: AbstractDict{SymbolicT, V}} <: AbstractDict{SymbolicT, V}
    dict::D

    function AtomicArrayDict(dict::AbstractDict{SymbolicT, V}) where {V}
        for k in keys(dict)
            validate_atomic_array_key(k)
        end
        new{V, typeof(dict)}(dict)
    end
end

AtomicArrayDict{V, D}(dict::AtomicArrayDict{V, D}) where {V, D} = copy(dict)
AtomicArrayDict{V, D}() where {V, D} = AtomicArrayDict(D())
AtomicArrayDict() = AtomicArrayDict(Dict{SymbolicT, SymbolicT}())
AtomicArrayDict(args::Pair...) = AtomicArrayDict(Dict(args...))
AtomicArrayDict{V, D}(args::Pair...) where {V, D} = AtomicArrayDict(Dict(args...))
AtomicArrayDict{V}(args...) where {V} = AtomicArrayDict(Dict{SymbolicT, V}(args...))
AtomicArrayDict{V, D}(args...) where {V, D} = AtomicArrayDict(D(args...))

struct IndexedArrayKeyError <: Exception
    k::SymbolicT
end

function Base.showerror(io::IO, err::IndexedArrayKeyError)
    print(io, """
    `AtomicArrayDict` treats symbolic arrays as atomic. It does not allow keys to be \
    indexed array symbolics. Got key $(err.k).
    """)
end

function validate_atomic_array_key(k::SymbolicT)
    split_indexed_var(k)[2] && throw(IndexedArrayKeyError(k))
end

Base.copy(dd::AtomicArrayDict) = AtomicArrayDict(copy(dd.dict))
function Base.empty(dd::AtomicArrayDict, ::Type{K}, ::Type{V}) where {K, V}
    AtomicArrayDict(empty(dd.dict, K, V))
end

Base.get(def::Base.Callable, dd::AtomicArrayDict, k) = def()
Base.get(def::Base.Callable, dd::AtomicArrayDict, k::SymbolicT) = get(def, dd.dict, k)
function Base.get(f::Base.Callable, dd::AtomicArrayDict, k::Union{Num, Arr, CallAndWrap})
    return get(f, dd, unwrap(k))
end
Base.get(dd::AtomicArrayDict, k, default) = get(Returns(default), dd, k)

Base.haskey(dd::AtomicArrayDict, k) = haskey(dd.dict, k)

Base.getindex(dd::AtomicArrayDict, k) = dd.dict[k]

function Base.setindex!(dd::AtomicArrayDict, v, k)
    k = unwrap(k)
    validate_atomic_array_key(unwrap(k))
    setindex!(dd.dict, v, k)
end

Base.isempty(dd::AtomicArrayDict) = isempty(dd.dict)
Base.length(dd::AtomicArrayDict) = length(dd.dict)
Base.iterate(dd::AtomicArrayDict, args...) = Base.iterate(dd.dict, args...)
Base.sizehint!(dd::AtomicArrayDict, n; kw...) = sizehint!(dd.dict, n; kw...)
Base.empty!(dd::AtomicArrayDict) = empty!(dd.dict)

Base.delete!(dd::AtomicArrayDict, k) = delete!(dd.dict, k)

struct AtomicArraySet{D <: AbstractDict{SymbolicT, Nothing}} <: AbstractSet{SymbolicT}
    dd::AtomicArrayDict{Nothing, D}

    function AtomicArraySet{D}(dd::AtomicArrayDict{Nothing, D}) where {D}
        new{D}(dd)
    end
end

AtomicArraySet() = AtomicArraySet{Dict{SymbolicT, Nothing}}()
AtomicArraySet{D}() where {D} = AtomicArraySet{D}(D())
AtomicArraySet{D}(x::D) where {D} = AtomicArraySet{D}(AtomicArrayDict(x))

Base.isempty(x::AtomicArraySet) = isempty(x.dd)
Base.length(x::AtomicArraySet) = length(x.dd)
Base.sizehint!(x::AtomicArraySet, n::Integer) = (sizehint!(x.dd, n); x)
Base.in(item, x::AtomicArraySet) = haskey(x.dd, item)
Base.push!(x::AtomicArraySet, item) = (x.dd[item] = nothing; x)
Base.delete!(x::AtomicArraySet, item) = (delete!(x.dd, item); x)
Base.empty(::AtomicArraySet{D}) where {D} = AtomicArraySet{D}()
Base.copy(x::AtomicArraySet{D}) where {D} = AtomicArraySet{D}(copy(x.dd))
Base.iterate(x::AtomicArraySet, args...) = iterate(keys(x.dd), args...)

function Base.filter!(f::F, x::AtomicArraySet) where {F}
    filter!(f ∘ first, x.dd)
    return x
end

"""
    $TYPEDSIGNATURES

Add `item` to `x`. If `item` is an indexed array, add the array instead.
"""
function push_as_atomic_array!(x::AtomicArraySet, item::SymbolicT)
    push!(x, split_indexed_var(item)[1])
end
