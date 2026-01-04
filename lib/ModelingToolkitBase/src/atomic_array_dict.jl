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
        return new{V, typeof(dict)}(dict)
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
    return print(
        io, """
        `AtomicArrayDict` treats symbolic arrays as atomic. It does not allow keys to be \
        indexed array symbolics. Got key $(err.k).
        """
    )
end

function validate_atomic_array_key(k::SymbolicT)
    return split_indexed_var(k)[2] && throw(IndexedArrayKeyError(k))
end

Base.copy(dd::AtomicArrayDict) = AtomicArrayDict(copy(dd.dict))
function Base.empty(dd::AtomicArrayDict, ::Type{K}, ::Type{V}) where {K, V}
    return AtomicArrayDict(empty(dd.dict, K, V))
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
    return setindex!(dd.dict, v, k)
end

Base.isempty(dd::AtomicArrayDict) = isempty(dd.dict)
Base.length(dd::AtomicArrayDict) = length(dd.dict)
Base.iterate(dd::AtomicArrayDict, args...) = Base.iterate(dd.dict, args...)
Base.sizehint!(dd::AtomicArrayDict, n; kw...) = sizehint!(dd.dict, n; kw...)
Base.empty!(dd::AtomicArrayDict) = empty!(dd.dict)

Base.delete!(dd::AtomicArrayDict, k) = delete!(dd.dict, k)

"""
    $TYPEDSIGNATURES

Convert the symbolic mapping `dict` to an `AtomicArrayDict`. If `dict` contains keys which
are elements of a symbolic array, the returned mappng will have a key for the array, and
a value which is a symbolic array where entries specified in `dict` are present and `default`
otherwise.
"""
function as_atomic_dict_with_defaults(dict::AbstractDict{SymbolicT, SymbolicT}, default::SymbolicT)
    dd = AtomicArrayDict(empty(dict))
    indexed_array_vals = empty(dict, SymbolicT, Array{SymbolicT})
    for (k, v) in dict
        arr, isarr = split_indexed_var(k)
        if isarr
            buffer = get!(() -> fill(default, size(arr)), indexed_array_vals, arr)
            si = get_stable_index(k)
            buffer[si] = v
        else
            dd[k] = v
        end
    end
    for (k, v) in indexed_array_vals
        if all(SU.isconst, v)
            dd[k] = BSImpl.Const{VartypeT}(unwrap_const.(v))
        else
            dd[k] = BSImpl.Const{VartypeT}(v)
        end
    end
    return dd
end

"""
    $TYPEDSIGNATURES

Modify an atomic array mapping `dd` to map `k` to `v`. If `k` is an indexed array symbolic,
update the array to have value `v` at the corresponding index. If the array is not a key,
create the key and set all other entries to `default`.
"""
function write_possibly_indexed_array!(dd::AtomicArrayDict{SymbolicT}, k::SymbolicT, v::SymbolicT, default::SymbolicT)
    arr, isarr = split_indexed_var(k)
    if isarr
        buffer::Array{SymbolicT} = if haskey(dd, arr)
            collect(dd[arr])
        else
            fill(default, size(arr))
        end
        idx = get_stable_index(k)
        buffer[idx] = v
        if all(SU.isconst, buffer)
            dd[arr] = BSImpl.Const{VartypeT}(unwrap_const.(buffer))
        else
            dd[arr] = BSImpl.Const{VartypeT}(buffer)
        end
    else
        dd[k] = v
    end
    return dd
end

"""
    $TYPEDSIGNATURES

Check if `dd` has the key `k`. If `k` is indexed, check if `dd` has the array as a key.
"""
function has_possibly_indexed_key(dd::AtomicArrayDict, k::SymbolicT)
    arr, _ = split_indexed_var(k)
    return haskey(dd, arr)
end

"""
    $TYPEDSIGNATURES

Equivalent to `get(dd, k, default)`. If `k` is an indexed array, then return
`dd[arr][idxs...]` for the corresponding array `arr` and indices, or `default`
if `arr` does not exist.
"""
function get_possibly_indexed(dd::AtomicArrayDict, k::SymbolicT, default)
    arr, isarr = split_indexed_var(k)
    res = get(dd, arr, default)
    isarr || return res
    res === default && return default
    idx = get_stable_index(k)
    return res[idx]
end

struct AtomicArraySet{D <: AbstractDict{SymbolicT, Nothing}} <: AbstractSet{SymbolicT}
    dd::AtomicArrayDict{Nothing, D}

    function AtomicArraySet{D}(dd::AtomicArrayDict{Nothing, D}) where {D}
        return new{D}(dd)
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
    return push!(x, split_indexed_var(item)[1])
end

"""
    $METHODLIST

Convert an array of possibly scalarized variables into an `AtomicArraySet`.
"""
as_atomic_array_set(vars::Vector{SymbolicT}) = as_atomic_array_set(Dict{SymbolicT, Nothing}, vars)
function as_atomic_array_set(::Type{D}, vars::Vector{SymbolicT}) where {D}
    set = AtomicArraySet{D}()
    for v in vars
        push_as_atomic_array!(set, v)
    end
    return set
end

function contains_possibly_indexed_element(x::AtomicArraySet, k::SymbolicT)
    return has_possibly_indexed_key(x.dd, k)
end

"""
    $TYPEDEF

A wrapper around `AtomicArrayDict{SymbolicT, D}` which makes it more amenable to
substitution. This wrapper allows indexing with indexed arrays, invoking
[`get_possibly_indexed`](@ref), [`has_possibly_indexed_key`](@ref) and
[`write_possibly_indexed_array!`](@ref) as appropriate. More significantly,
when `Base.get` (and `Base.getindex`) is called with an indexed array for which
the corresponding array is present in the wrapped `dict` but the value is `default`,
it will return the indexed array. For example, if the wrapped `dict` has
`k => [default, other_var]` as a key-value pair and this wrapper is accessed at
`k[1]`, it will return `k[1]` instead of `default`. This is useful since `default`
is used to represent missing values. Substituting something like `sin(k[1])` will
then not attempt to perform `sin(default)` but instead return `sin(k[1])`.
"""
struct AtomicArrayDictSubstitutionWrapper{D} <: AbstractDict{SymbolicT, SymbolicT}
    dict::AtomicArrayDict{SymbolicT, D}
    default::SymbolicT
end

function AtomicArrayDictSubstitutionWrapper(d::AtomicArrayDict{SymbolicT, D}) where {D}
    return AtomicArrayDictSubstitutionWrapper{D}(d, COMMON_NOTHING)
end

const AADSubWrapper{D} = AtomicArrayDictSubstitutionWrapper{D}

Base.get(def::Base.Callable, dd::AADSubWrapper, k) = def()
function Base.get(def::Base.Callable, dd::AADSubWrapper, k::SymbolicT)
    arr, isarr = split_indexed_var(k)
    res = get_possibly_indexed(dd.dict, k, dd.default)
    if res === dd.default
        isarr && haskey(dd.dict, arr) && return k
        return def()
    end
    return res
end
function Base.get(f::Base.Callable, dd::AADSubWrapper, k::Union{Num, Arr, CallAndWrap})
    return get(f, dd, unwrap(k))
end
Base.get(dd::AADSubWrapper, k, default) = get(Returns(default), dd, k)

function Base.haskey(dd::AADSubWrapper, k::SymbolicT)
    return has_possibly_indexed_key(dd, k)
end
function Base.haskey(dd::AADSubWrapper, k::Union{Num, Arr, CallAndWrap})
    return haskey(dd, unwrap(k))
end

function Base.getindex(dd::AADSubWrapper, k)
    res = get(dd, k, dd.default)
    res === dd.default && throw(KeyError(k))
    return res
end

function Base.setindex!(dd::AADSubWrapper, v, k)
    return write_possibly_indexed_array!(dd.dict, k, v, dd.default)
end

Base.isempty(dd::AADSubWrapper) = isempty(dd.dict)
Base.length(dd::AADSubWrapper) = length(dd.dict)
Base.iterate(dd::AADSubWrapper, args...) = Base.iterate(dd.dict, args...)
Base.sizehint!(dd::AADSubWrapper, n; kw...) = sizehint!(dd.dict, n; kw...)
Base.empty!(dd::AADSubWrapper) = empty!(dd.dict)
