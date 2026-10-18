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

struct RecordFieldKeyError <: Exception
    k::SymbolicT
end

function Base.showerror(io::IO, err::RecordFieldKeyError)
    return print(
        io, """
        `AtomicArrayDict` treats symbolic records as atomic. It does not allow keys to be \
        projections of a record symbolic. Got key $(err.k).
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
    while true
        Moshi.Match.@match x begin
            BSImpl.Term(; f, args) && if f isa Symbolics.SymbolicGetproperty || f === getindex end => begin
                x = args[1]::SymbolicT
            end
            _ => break
        end
    end
    return x, Symbolics.issymstruct(x)
end

function validate_atomic_array_key(k::SymbolicT)
    # Record projections are not rejected here. Initialization is keyed on leaves -
    # `Initial` exists for `p.x` and not for `p`, a record having no shape in the numeric
    # initials buffer - and the initialization system it builds is itself a `System`, whose
    # stored maps are of this type. What a user hands in is checked on the way in, by
    # [`validate_given_atomic_array_key`](@ref).
    return split_indexed_var(k)[2] && throw(IndexedArrayKeyError(k))
end

"""
    $TYPEDSIGNATURES

Check a key of a map a user supplied, which is stricter than
[`validate_atomic_array_key`](@ref): a projection of a record is rejected the way an array
element is, because the value belongs in the whole record, which is what gets stored.
Setting one field is spelled `write_possibly_indexed_array!`, exactly as for an array.
"""
function validate_given_atomic_array_key(k::SymbolicT)
    root, kind = split_composite_var(k)
    kind === COMPOSITE_ARRAY && throw(IndexedArrayKeyError(k))
    kind === COMPOSITE_RECORD && !isequal(root, k) && throw(RecordFieldKeyError(k))
    return false
end

"""
    $TYPEDSIGNATURES

The leaf-wise value buffer stored for a record key.
"""
record_buffer(stored::SymbolicT) = collect(Iterators.map(Base.Fix1(getindex, stored), SU.stable_eachindex(stored)))

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
    # `root => [(index => value), ...]` for every entry written through an array element.
    indexed_array_writes = empty(dict, SymbolicT, Vector{Pair{Any, SymbolicT}})
    record_vals = empty(dict, SymbolicT, Vector{SymbolicT})
    # `(root, node, value, leaf indices)` for every entry stored under a record.
    record_writes = Tuple{SymbolicT, SymbolicT, SymbolicT, Vector{Int}}[]
    for (k, v) in dict
        root, kind = split_composite_var(k)
        if kind === COMPOSITE_RECORD
            push!(record_writes, (root, k, v, record_node_indices(root, k)))
        elseif kind === COMPOSITE_ARRAY
            # Collected rather than applied here: the buffer they are applied to depends
            # on whether `dict` also writes the whole array, which may come later in
            # iteration order.
            push!(get!(() -> Pair{Any, SymbolicT}[], indexed_array_writes, root),
                get_stable_index(k) => v)
        else
            dd[k] = v
        end
    end
    # Broadest writes first, as for records below: an entry for one element overrides the
    # whole-array entry it belongs to, rather than replacing it wholesale. This matches
    # `write_possibly_indexed_array!`, which seeds from the existing entry.
    for (root, writes) in indexed_array_writes
        buffer::Array{SymbolicT} = if haskey(dd, root)
            collect(dd[root])
        else
            fill(default, size(root))
        end
        isempty(buffer) && continue
        for (si, v) in writes
            buffer[si] = v
        end
        if all(SU.isconst, buffer)
            __unsafe_aad_setindex!(dd, BSImpl.Const{VartypeT}(unwrap_const.(buffer)), root)
        else
            __unsafe_aad_setindex!(dd, BSImpl.Const{VartypeT}(buffer), root)
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
    for (root, buffer) in record_vals
        __unsafe_aad_setindex!(dd, record_value_from_buffer(root, buffer), root)
    end
    return dd
end

"""
    $TYPEDSIGNATURES

The value to store for record `root` given `buffer`, its leaves in
[`record_leaves`](@ref) order. Unfilled leaves are still `default`, so the result is a
symbolic literal until every one of them has been written - the record counterpart of an
array entry that still holds the sentinel at some index.
"""
function record_value_from_buffer(root::SymbolicT, buffer::Vector{SymbolicT})
    vals = map(x -> SU.isconst(x) ? unwrap_const(x) : x, buffer)
    return BSImpl.Const{VartypeT}(record_node_value(root, root, vals))
end

"""
    $TYPEDSIGNATURES

The value to store for record `root` given one value per [`record_lowering_nodes`](@ref).

Those nodes are not the leaves: a numeric array field is one node covering several leaves.
Its value is spread back over them, and a field left unfilled marks every leaf beneath it
unfilled too.
"""
function record_value_from_nodes(
        root::SymbolicT, nodes::Vector{SymbolicT}, nodevals::Vector{SymbolicT}
    )
    leaves = record_leaves(root)
    buffer = SymbolicT[COMMON_NOTHING for _ in eachindex(leaves)]
    for (node, nv) in zip(nodes, nodevals)
        idxs = record_node_indices(root, node)
        if length(idxs) == 1
            buffer[only(idxs)] = nv
        else
            for (j, i) in enumerate(idxs)
                buffer[i] = nv === COMMON_NOTHING ? COMMON_NOTHING : nv[j]
            end
        end
    end
    return record_value_from_buffer(root, buffer)
end

"""
    $TYPEDSIGNATURES

Modify an atomic array mapping `dd` to map `k` to `v`. If `k` is an indexed array symbolic,
update the array to have value `v` at the corresponding index. If the array is not a key,
create the key and set all other entries to `default`. If `k` is an array symbolic and `v`
is scalar-shaped, `v` is broadcast to `size(k)`.
"""
function write_possibly_indexed_array!(dd::AtomicArrayDict{SymbolicT}, k::SymbolicT, v::SymbolicT, default::SymbolicT)
    rroot, isrec = record_key_root(k)
    # A projection is folded into the record it belongs to, the way an array element is
    # folded into its array - unless the map is already keyed on that projection, which is
    # how initialization's maps are keyed. Writing the root there would leave the value in
    # two places, and the stale leaf is the one a lookup finds.
    if isrec && !isequal(rroot, k) && !haskey(dd, k)
        return write_record_field!(dd, rroot, k, v, default)
    end
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
        shk = SU.shape(k)
        dd[k] = if SU.is_array_shape(shk) && !(shk isa SU.Unknown) &&
                !SU.is_array_shape(SU.shape(v))
            if SU.isconst(v)
                BSImpl.Const{VartypeT}(fill(unwrap_const(v), size(k)))
            else
                BSImpl.Const{VartypeT}(fill(v, size(k)))
            end
        else
            v
        end
    end
    return dd
end

"""
    $TYPEDSIGNATURES

Fill `v` into field `k` of the record stored under `root`, leaving its other fields as they
were - the record counterpart of writing one element of a stored array.

The record is read back into its leaves, the leaves this field covers are written, and the
whole is reassembled and stored again. Leaves never written stay `default`, so a record
that has only been partly specified reassembles into a symbolic literal rather than a value
of the struct type, which has no way to say "not yet known".
"""
function write_record_field!(
        dd::AtomicArrayDict{SymbolicT}, root::SymbolicT, k::SymbolicT, v::SymbolicT,
        default::SymbolicT
    )
    leaves = record_leaves(root)
    buffer::Vector{SymbolicT} = if haskey(dd, root)
        # Read the stored record back into its leaves by projecting each one out of it,
        # the same way a stored array is read back with `collect`.
        stored = dd[root]
        SymbolicT[record_leaf_entry(leaf, root, stored) for leaf in leaves]
    else
        fill(default, length(leaves))
    end
    isempty(buffer) && return dd
    for i in record_node_indices(root, k)
        buffer[i] = record_leaf_entry(leaves[i], k, v)
    end
    __unsafe_aad_setindex!(dd, record_value_from_buffer(root, buffer), root)
    return dd
end

"""
    $TYPEDSIGNATURES

Merge `k => v` into `dd` the way [`left_merge!`](@ref) would, for a `k` that projects out of
a record. The record is stored whole, so this is a partial write of the fields `k` covers,
skipping the ones `dd` already has a value for - the same "left wins, except where left
doesn't know" rule, applied through the record rather than to a key of its own.
"""
function merge_record_field!(dd::AtomicArrayDict{SymbolicT}, k::SymbolicT, v::SymbolicT)
    if SU.is_array_shape(SU.shape(k))
        for i in SU.stable_eachindex(k)
            ki = k[i]::SymbolicT
            vi = v[i]
            vi === COMMON_NOTHING && continue
            get_possibly_indexed(dd, ki, COMMON_NOTHING) === COMMON_NOTHING || continue
            write_possibly_indexed_array!(dd, ki, vi, COMMON_NOTHING)
        end
    elseif get_possibly_indexed(dd, k, COMMON_NOTHING) === COMMON_NOTHING
        write_possibly_indexed_array!(dd, k, v, COMMON_NOTHING)
    end
    return dd
end

"""
    $TYPEDSIGNATURES

The nodes a record key expands to for initialization: the deepest projections of `root`
whose `Initial` exists. That is a leaf, except where a field is a numeric array, which
stays whole - an array has an `Initial` of its own, so expanding past it would ask for a
parameter that was never created.
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

Replace every record key of `op` with the nodes `Initial` is defined for, each holding
that node's projection of
the record's value. `op` is modified, since callers of this go on to write their own
entries into the same map and expect to see them afterwards.

Records are stored atomically, but initialization speaks in leaves: `Initial` exists for
`p.x` and not for `p`, since a record has no shape to occupy in the flat numeric initials
buffer - unlike an array, which gets an `Initial` of its own and so needs no translation.
This converts the one vocabulary into the other at the door of initialization, leaving
storage untouched; an unfilled field projects to the same sentinel a partially written
array element would hold.
"""
function expand_record_keys!(op::AbstractDict)
    records = SymbolicT[k for k in keys(op) if Symbolics.issymstruct(k)]
    isempty(records) && return op
    for k in records
        v = op[k]
        delete!(op, k)
        for node in record_lowering_nodes(k)
            entry = record_leaf_entry(node, k, v)
            # A field the record says nothing about contributes no entry, exactly as a
            # scalar variable nobody gave a value for is simply absent from the map. An
            # entry holding "unfilled" would instead read as a value, and would shadow the
            # one initialization goes on to determine.
            is_unfilled_record_entry(entry) && continue
            __unsafe_aad_setindex!(op, entry, node)
        end
    end
    return op
end

"""
    $TYPEDSIGNATURES

Whether a value projected out of a record says nothing about that field. The sentinels
standing for "not given" are the cases; `CommonSentinel` is among them because
`build_operating_point` substitutes it for an explicit `nothing`.
"""
function is_unfilled_record_entry(v::SymbolicT)
    v === COMMON_NOTHING && return true
    v === COMMON_MISSING && return true
    v === COMMON_SENTINEL && return true
    return SU.isconst(v) && unwrap_const(v) isa CommonSentinel
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
projection is found under the record it belongs to, which is stored atomically and is not
necessarily its immediate parent.
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

When `k` is a scalar `getindex` term (possibly under unary [`Operator`](@ref)s)
with non-constant indices and `arr` is present:
- if `allow_symbolic_indices` is `false` (default), return `default` so
  membership-style callers treat the entry as missing until indices are concrete;
- if `allow_symbolic_indices` is `true`, return the symbolic `getindex` of the
  substituted array value so a fixpoint substituter can resolve the indices next.

Slice / `ArrayOp` and other non-`getindex` forms always go through
[`get_stable_index`](@ref), which throws.
"""
function get_possibly_indexed(
        dd::AtomicArrayDict, k::SymbolicT, default; allow_symbolic_indices::Bool = false
    )
    arr, isarr = split_indexed_var(k)
    res = get(dd, arr, default)
    if !isarr
        res === default || return res
        return record_ancestor_value(dd, k, default)
    end
    res === default && return record_ancestor_value(dd, k, default)
    if !has_const_int_indices(k) && is_getindex_indexed(k)
        allow_symbolic_indices || return default
        return index_substituted_array(res, k)
    end
    idx = get_stable_index(k)
    return res[idx]
end

"""
    $TYPEDSIGNATURES

The value `k` projects out of whichever ancestor of it `dd` holds, or `default` if it holds
none. Records are stored atomically, so a field is found under the record it belongs to.
"""
function record_ancestor_value(dd::AtomicArrayDict, k::SymbolicT, default)
    anc = nearest_present_ancestor(dd, k)
    anc === nothing && return default
    return record_leaf_entry(k, anc, dd[anc])
end



"""
    $TYPEDSIGNATURES

Return `true` if `k` is a scalar `getindex` term, possibly under unary
[`Operator`](@ref) wrappers. Slice / [`ArrayOp`](@ref) forms return `false`.
"""
function is_getindex_indexed(k::SymbolicT)
    return Moshi.Match.@match k begin
        BSImpl.Term(; f, args) && if f === getindex end => true
        BSImpl.Term(; f, args) && if f isa Operator && length(args) == 1 end => begin
            return is_getindex_indexed(args[1]::SymbolicT)
        end
        _ => return false
    end
end

"""
    $TYPEDSIGNATURES

Build `getindex(res, idxs...)` using the (possibly symbolic) indices of the
indexed term `k`. Recurses through wrapping [`Operator`](@ref)s.
"""
function index_substituted_array(res::SymbolicT, k::SymbolicT)
    return Moshi.Match.@match k begin
        BSImpl.Term(; f, args) && if f === getindex end => begin
            return res[args[2:end]...]::SymbolicT
        end
        BSImpl.Term(; f, args) && if f isa Operator && length(args) == 1 end => begin
            # `res` is already the value of `f(arr)`; do not re-apply `f`
            # (mirrors `_get_stable_index`, which discards the operator).
            return index_substituted_array(res, args[1]::SymbolicT)::SymbolicT
        end
        _ => throw(ArgumentError(lazy"Invalid variable $k for `index_substituted_array`."))
    end
end

"""
    $(TYPEDEF)

A set of symbolic variables, backed by an [`AtomicArrayDict`](@ref).

It holds the two kinds of composite at different granularities, which is deliberate:

  - A symbolic array is one entry. [`push_as_atomic_array!`](@ref) normalizes `x[1]` to
    `x`, matching the dictionary, where `x` is a legal key and `x[1]` is not.
  - A record is held at whatever granularity it is pushed at, and in practice that is its
    leaves. A record's unknowns genuinely are its leaves - `p.x` is solved for, while `p`
    is not a quantity with a value of its own - and the sets backed by this type are built
    from maps that [`expand_record_keys!`](@ref) has already translated into leaves.

Membership is therefore exact: a set given `p` does not report `p.x` as a member, and one
given the leaves does not report `p`. Nothing here stops a caller from putting both
spellings into one set. That the current callers do not follows from the order their maps
are built, not from an invariant this type enforces, so a new caller that pushes a record
whole alongside leaf-wise ones would get a set whose membership answers disagree.
"""
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
# Deliberately unchecked: this is where a record goes in at whatever granularity the
# caller has, which the docstring above describes. Array keys are normalized by
# `push_as_atomic_array!` instead, before they reach here.
Base.push!(x::AtomicArraySet, item) = (__unsafe_aad_setindex!(x.dd, nothing, unwrap(item)); x)
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

Add `item` to `x`. If `item` is an indexed array, add the array instead. A record
projection is added as-is, at leaf granularity - see [`AtomicArraySet`](@ref) for why the
two composites are treated differently.
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

Limitation: when the stored array value uses `COMMON_NOTHING` holes (e.g.
`Union{Nothing,Float64}[1.0, nothing, 3.0]`), looking up a *symbolic* index such
as `a[i]` is not supported — constructing the intermediate symbolic `getindex`
of that heterogeneous array currently errors. Concrete integer indices still
honour the hole contract above. Prefer leaving those slots unset in the map
until indices are concrete, or fill holes before substitution.
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
    res = get_possibly_indexed(dd.dict, k, dd.default; allow_symbolic_indices = true)
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
