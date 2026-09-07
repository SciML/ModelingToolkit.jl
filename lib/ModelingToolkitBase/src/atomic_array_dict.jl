"""
    $(TYPEDEF)

Wrapper over an `AbstractDict{SymbolicT, V} where {V}` which disallows keys that are
indexed array symbolics. Specifically, if `@variables x[1:4]` exists, then `x` can be
a key but `x[1]` cannot.
"""
struct AtomicArrayDict{V, D <: AbstractDict{SymbolicT, V}} <: AbstractDict{SymbolicT, V}
    dict::D

    function AtomicArrayDict(dict::AbstractDict{SymbolicT, V}; __check = true) where {V}
        if __check
            for k in keys(dict)
                validate_atomic_array_key(k)
            end
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

"""
    $TYPEDEF

How a symbolic key relates to the composite it is stored under. `COMPOSITE_NONE` means the
key is atomic in its own right.
"""
@enum CompositeKind COMPOSITE_NONE COMPOSITE_ARRAY COMPOSITE_RECORD

"""
    $TYPEDSIGNATURES

Given a symbolic `k`, return the composite it is stored under and which kind it is. Struct
symbolics are stored under the whole record, array elements under the whole array. A key
that is its own root is `COMPOSITE_NONE`, except for a whole record, which is
`COMPOSITE_RECORD` because its value is laid out leaf-wise.
"""
function split_composite_var(k::SymbolicT)::Tuple{SymbolicT, CompositeKind}
    root, isrec = record_key_root(k)
    isrec && return root, COMPOSITE_RECORD
    arr, isarr = split_indexed_var(k)
    isarr && return arr, COMPOSITE_ARRAY
    return k, COMPOSITE_NONE
end

"""
    $TYPEDSIGNATURES

Peel field accesses and indices off `k` and return the result, along with whether it is a
struct symbolic. Unlike [`split_field_access`](@ref) this does *not* commute with
operators: `Initial(p.y)` is its own key, not a projection of `Initial(p)`. Operators are
applied to leaves, not to records, since a record symtype does not flatten into a numeric
parameter buffer.
"""
function record_key_root(k::SymbolicT)::Tuple{SymbolicT, Bool}
    x = k
    while iscall(x)
        f = operation(x)
        f isa Symbolics.SymbolicGetproperty || f === getindex || break
        x = arguments(x)[1]::SymbolicT
    end
    return x, Symbolics.issymstruct(x)
end

function validate_atomic_array_key(k::SymbolicT)
    return split_indexed_var(k)[2] && throw(IndexedArrayKeyError(k))
end

"""
    $TYPEDSIGNATURES

The leaf-wise value buffer stored for a record key.
"""
record_buffer(stored::SymbolicT) = Vector{SymbolicT}(collect(stored))

Base.copy(dd::AtomicArrayDict) = AtomicArrayDict(copy(dd.dict); __check = false)
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
    validate_atomic_array_key(k)
    return __unsafe_aad_setindex!(dd, v, k)
end

function __unsafe_aad_setindex!(dd::AtomicArrayDict, v, k::SymbolicT)
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
    record_vals = empty(dict, SymbolicT, Vector{SymbolicT})
    # `(root, node, value, leaf indices)` for every entry stored under a record.
    record_writes = Tuple{SymbolicT, SymbolicT, SymbolicT, Vector{Int}}[]
    for (k, v) in dict
        root, kind = split_composite_var(k)
        if kind === COMPOSITE_RECORD
            push!(record_writes, (root, k, v, record_node_indices(root, k)))
        elseif kind === COMPOSITE_ARRAY
            buffer = get!(() -> fill(default, size(root)), indexed_array_vals, root)
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
    # Broadest writes first, so that a more specific entry overrides the record it
    # belongs to rather than the other way round.
    sort!(record_writes; by = length ∘ last, rev = true)
    for (root, node, v, idxs) in record_writes
        leaves = record_leaves(root)
        buffer = get!(() -> fill(default, length(leaves)), record_vals, root)
        for i in idxs
            buffer[i] = record_leaf_entry(leaves[i], node, v)
        end
    end
    for (k, v) in record_vals
        __unsafe_aad_setindex!(dd, BSImpl.Const{VartypeT}(v), k)
    end
    # Records are the right granularity for resolving what was asked for, but nothing
    # downstream expects a record key, so lower before returning.
    lower_record_entries!(dd, default)
    return dd
end

"""
    $TYPEDSIGNATURES

Modify an atomic array mapping `dd` to map `k` to `v`. If `k` is an indexed array symbolic,
update the array to have value `v` at the corresponding index. If the array is not a key,
create the key and set all other entries to `default`.
"""
function write_possibly_indexed_array!(dd::AtomicArrayDict{SymbolicT}, k::SymbolicT, v::SymbolicT, default::SymbolicT)
    root, isarr = split_indexed_var(k)
    if isarr
        buffer::Array{SymbolicT} = if haskey(dd, root)
            collect(dd[root])
        else
            fill(default, size(root))
        end
        isempty(buffer) && return dd
        idx = get_stable_index(k)
        buffer[idx] = v
        if all(SU.isconst, buffer)
            dd[root] = BSImpl.Const{VartypeT}(unwrap_const.(buffer))
        else
            dd[root] = BSImpl.Const{VartypeT}(buffer)
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
    haskey(dd, k) && return true
    haskey(dd, split_indexed_var(k)[1]) && return true
    return nearest_present_ancestor(dd, k) !== nothing
end

"""
    $TYPEDSIGNATURES

The innermost proper ancestor of `k` that is a key of `dd`, or `nothing`. A struct
projection is stored under whichever ancestor `lower_record_entries!` stopped at, which is
not necessarily its immediate parent.
"""
function nearest_present_ancestor(dd::AtomicArrayDict, k::SymbolicT)
    x = k
    while iscall(x)
        f = operation(x)
        f isa Symbolics.SymbolicGetproperty || f === getindex || break
        x = arguments(x)[1]::SymbolicT
        haskey(dd, x) && return x
    end
    return nothing
end

"""
    $TYPEDSIGNATURES

Equivalent to `get(dd, k, default)`. If `k` is an indexed array, then return
`dd[arr][idxs...]` for the corresponding array `arr` and indices, or `default`
if `arr` does not exist.
"""
function get_possibly_indexed(dd::AtomicArrayDict, k::SymbolicT, default)
    haskey(dd, k) && return dd[k]
    arr, isarr = split_indexed_var(k)
    if isarr
        res = get(dd, arr, default)
        res === default || return res[get_stable_index(k)]
    end
    anc = nearest_present_ancestor(dd, k)
    anc === nothing && return default
    return record_leaf_entry(k, anc, dd[anc])
end

"""
    $TYPEDSIGNATURES

The nodes a record key is lowered to: the deepest projections of `root` whose symtype
flattens into a numeric buffer. This is a leaf, except where a field is a numeric array,
which stays whole because ordinary array atomicity already covers it. It is the same
granularity `Initial` uses, which is what the initialization system pairs `op` keys with.
"""
function record_lowering_nodes(root::SymbolicT)::Vector{SymbolicT}
    nodes = SymbolicT[]
    _record_lowering_nodes!(nodes, root)
    return nodes
end

function _record_lowering_nodes!(nodes::Vector{SymbolicT}, node::SymbolicT)
    T = SU.symtype(node)
    if !Symbolics.issymstruct(node)
        if SU.is_array_shape(SU.shape(node)) && !(T <: AbstractArray{<:Number})
            for idx in SU.stable_eachindex(node)
                _record_lowering_nodes!(nodes, node[idx]::SymbolicT)
            end
            return nodes
        end
        push!(nodes, node)
        return nodes
    end
    for fname in fieldnames(T::DataType)
        field = Symbolics.SymbolicGetproperty{T, fname}()(node)::SymbolicT
        _record_lowering_nodes!(nodes, field)
    end
    return nodes
end

"""
    $TYPEDSIGNATURES

Replace every record key of `dd` with its [`record_lowering_nodes`](@ref). Records are the
right granularity for resolving what the user asked for, but everything downstream of that
is keyed per numeric buffer entry. Entries still holding `default` are dropped.
"""
function lower_record_entries!(dd::AtomicArrayDict{SymbolicT}, default::SymbolicT)
    for k in collect(keys(dd))
        Symbolics.issymstruct(k) || continue
        leaves = record_leaves(k)
        buffer = record_buffer(dd[k])
        lowered = Pair{SymbolicT, SymbolicT}[]
        for node in record_lowering_nodes(k)
            idxs = record_node_indices(k, node)
            any(i -> buffer[i] === default, idxs) && continue
            if length(idxs) == 1 && isequal(node, leaves[only(idxs)])
                push!(lowered, node => buffer[only(idxs)])
                continue
            end
            # An intermediate node is only known once every leaf below it is.
            all(SU.isconst, view(buffer, idxs)) || continue
            vals = map(x -> SU.isconst(x) ? unwrap_const(x) : x, buffer)
            push!(lowered, node => BSImpl.Const{VartypeT}(record_node_value(k, node, vals)))
        end
        delete!(dd, k)
        for (node, v) in lowered
            __unsafe_aad_setindex!(dd, v, node)
        end
    end
    return dd
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

Add `item` to `x`. If `item` is an indexed array, add the array instead. Struct projections
are added as-is: these sets track which variables exist, which is leaf granularity.
"""
function push_as_atomic_array!(x::AtomicArraySet, item::SymbolicT)
    return push!(x, split_indexed_var(item)[1])
end

"""
    $METHODLIST

Convert an array of possibly scalarized variables into an `AtomicArraySet`.
"""
as_atomic_array_set(vars::Vector{SymbolicT}) = as_atomic_array_set(Dict{SymbolicT, Nothing}, vars)
as_atomic_array_set(vars::Vector) = as_atomic_array_set(Vector{SymbolicT}(map(unwrap, vars)))
as_atomic_array_set(vars::AtomicArraySet) = vars
function as_atomic_array_set(
        ::Type{D},
        vars::Union{AbstractVector{SymbolicT}, AbstractSet{SymbolicT}}
    ) where {D}
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
    res = get_possibly_indexed(dd.dict, k, dd.default)
    if res === dd.default
        arr, isarr = split_indexed_var(k)
        isarr && haskey(dd.dict, arr) && return k
        nearest_present_ancestor(dd.dict, k) === nothing || return k
        return def()
    end
    if SU.is_array_shape(SU.shape(k)) && any(
            Base.Fix2(===, dd.default) ∘ Base.Fix1(getindex, res),
            SU.stable_eachindex(res)
        )
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
