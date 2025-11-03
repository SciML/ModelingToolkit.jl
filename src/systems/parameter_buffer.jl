symconvert(::Type{Symbolics.Struct{T}}, x) where {T} = convert(T, x)
symconvert(::Type{T}, x::V) where {T, V} = convert(promote_type(T, V), x)
symconvert(::Type{Real}, x::Integer) = convert(Float16, x)
symconvert(::Type{V}, x) where {V <: AbstractArray} = convert(V, symconvert.(eltype(V), x))

struct MTKParameters{T, I, D, C, N, H}
    tunable::T
    initials::I
    discrete::D
    constant::C
    nonnumeric::N
    caches::H
end

"""
    function MTKParameters(sys::AbstractSystem, p, u0 = Dict(); t0 = nothing)

Create an `MTKParameters` object for the system `sys`. `p` (`u0`) are symbolic maps from
parameters (unknowns) to their values. The values can also be symbolic expressions, which
are evaluated given the values of other parameters/unknowns. `u0` is only required if
the values of parameters depend on the unknowns. `t0` is the initial time, for time-
dependent systems. It is only required if the symbolic expressions also use the independent
variable of the system.

This requires that `complete` has been called on the system (usually via
`mtkcompile` or `@mtkcompile`) and the keyword `split = true` was passed (which is
the default behavior).
"""
function MTKParameters(
        sys::AbstractSystem, op; tofloat = false,
        t0 = nothing, substitution_limit = 1000, floatT = nothing,
        p_constructor = identity, fast_path = false)
    ic = if has_index_cache(sys) && get_index_cache(sys) !== nothing
        get_index_cache(sys)
    else
        error("Cannot create MTKParameters if system does not have index_cache")
    end
    all_ps = Set(unwrap.(parameters(sys; initial_parameters = true)))
    union!(all_ps, default_toterm.(unwrap.(parameters(sys; initial_parameters = true))))

    dvs = unknowns(sys)
    ps = parameters(sys; initial_parameters = true)
    op = to_varmap(op, ps)
    symbols_to_symbolics!(sys, op)
    defs = add_toterms(recursive_unwrap(defaults(sys)))

    is_time_dependent(sys) && add_observed!(sys, op)
    add_parameter_dependencies!(sys, op)

    u0map = anydict()
    pmap = anydict()
    if fast_path
        missing_pars = missingvars(op, ps)
    else
        _, missing_pars = build_operating_point!(sys, op,
            u0map, pmap, defs, dvs, ps)
    end
    if t0 !== nothing
        op[get_iv(sys)] = t0
    end

    if floatT === nothing
        floatT = float(float_type_from_varmap(op))
    end

    isempty(missing_pars) || throw(MissingParametersError(collect(missing_pars)))
    evaluate_varmap!(op, ps; limit = substitution_limit)

    p = op
    filter!(p) do kvp
        kvp[1] in all_ps
    end

    tunable_buffer = Vector{ic.tunable_buffer_size.type}(
        undef, ic.tunable_buffer_size.length)
    initials_buffer = Vector{ic.initials_buffer_size.type}(
        undef, ic.initials_buffer_size.length)
    disc_buffer = Tuple(BlockedArray(
                            Vector{subbuffer_sizes[1].type}(
                                undef, sum(x -> x.length, subbuffer_sizes)),
                            map(x -> x.length, subbuffer_sizes))
    for subbuffer_sizes in ic.discrete_buffer_sizes)
    const_buffer = Tuple(Vector{temp.type}(undef, temp.length)
    for temp in ic.constant_buffer_sizes)
    nonnumeric_buffer = Tuple(Vector{temp.type}(undef, temp.length)
    for temp in ic.nonnumeric_buffer_sizes)
    function set_value(sym, val)
        done = true
        if haskey(ic.tunable_idx, sym)
            idx = ic.tunable_idx[sym]
            tunable_buffer[idx] = val
        elseif haskey(ic.initials_idx, sym)
            idx = ic.initials_idx[sym]
            initials_buffer[idx] = val
        elseif haskey(ic.discrete_idx, sym)
            idx = ic.discrete_idx[sym]
            disc_buffer[idx.buffer_idx][idx.idx_in_buffer] = val
        elseif haskey(ic.constant_idx, sym)
            i, j = ic.constant_idx[sym]
            const_buffer[i][j] = val
        elseif haskey(ic.nonnumeric_idx, sym)
            i, j = ic.nonnumeric_idx[sym]
            nonnumeric_buffer[i][j] = val
        elseif !isequal(default_toterm(sym), sym)
            done = set_value(default_toterm(sym), val)
        else
            done = false
        end
        return done
    end
    for (sym, val) in p
        sym = unwrap(sym)
        val = unwrap(val)
        ctype = symtype(sym)
        if symbolic_type(val) !== NotSymbolic()
            error("Could not evaluate value of parameter $sym. Missing values for variables in expression $val.")
        end
        if ctype <: FnType
            ctype = fntype_to_function_type(ctype)
        end
        if ctype == Real && floatT !== nothing
            ctype = floatT
        end
        val = symconvert(ctype, val)
        done = set_value(sym, val)
        if !done && Symbolics.isarraysymbolic(sym)
            if Symbolics.shape(sym) === Symbolics.Unknown()
                for i in eachindex(val)
                    set_value(sym[i], val[i])
                end
            else
                if size(sym) != size(val)
                    error("Got value of size $(size(val)) for parameter $sym of size $(size(sym))")
                end
                set_value.(collect(sym), val)
            end
        end
    end
    tunable_buffer = narrow_buffer_type(tunable_buffer; p_constructor)
    if isempty(tunable_buffer)
        tunable_buffer = SizedVector{0, Float64}()
    end
    initials_buffer = narrow_buffer_type(initials_buffer; p_constructor)
    if isempty(initials_buffer)
        initials_buffer = SizedVector{0, Float64}()
    end
    disc_buffer = narrow_buffer_type.(disc_buffer; p_constructor)
    const_buffer = narrow_buffer_type.(const_buffer; p_constructor)
    # Don't narrow nonnumeric types
    if !isempty(nonnumeric_buffer)
        nonnumeric_buffer = map(p_constructor, nonnumeric_buffer)
    end

    mtkps = MTKParameters{
        typeof(tunable_buffer), typeof(initials_buffer), typeof(disc_buffer),
        typeof(const_buffer), typeof(nonnumeric_buffer), typeof(())}(tunable_buffer,
        initials_buffer, disc_buffer, const_buffer, nonnumeric_buffer, ())
    return mtkps
end

function rebuild_with_caches(p::MTKParameters, cache_templates::BufferTemplate...)
    buffers = map(cache_templates) do template
        Vector{template.type}(undef, template.length)
    end
    @set p.caches = buffers
end

function narrow_buffer_type(buffer::AbstractArray; p_constructor = identity)
    type = Union{}
    for x in buffer
        type = promote_type(type, typeof(x))
    end
    return p_constructor(type.(buffer))
end

function narrow_buffer_type(
        buffer::AbstractArray{<:AbstractArray}; p_constructor = identity)
    type = Union{}
    for arr in buffer
        for x in arr
            type = promote_type(type, typeof(x))
        end
    end
    buffer = map(buffer) do buf
        p_constructor(type.(buf))
    end
    return p_constructor(buffer)
end

function narrow_buffer_type(buffer::BlockedArray; p_constructor = identity)
    if eltype(buffer) <: AbstractArray
        buffer = narrow_buffer_type.(buffer; p_constructor)
    end
    type = Union{}
    for x in buffer
        type = promote_type(type, typeof(x))
    end
    tmp = p_constructor(type.(buffer))
    blocks = ntuple(Val(ndims(buffer))) do i
        bsizes = blocksizes(buffer, i)
        p_constructor(Int.(bsizes))
    end
    return BlockedArray(tmp, blocks...)
end

function buffer_to_arraypartition(buf)
    return ArrayPartition(ntuple(i -> _buffer_to_arrp_helper(buf[i]), Val(length(buf))))
end

_buffer_to_arrp_helper(v::T) where {T} = _buffer_to_arrp_helper(eltype(T), v)
_buffer_to_arrp_helper(::Type{<:AbstractArray}, v) = buffer_to_arraypartition(v)
_buffer_to_arrp_helper(::Any, v) = v

function _split_helper(buf_v::T, recurse, raw, idx) where {T}
    _split_helper(eltype(T), buf_v, recurse, raw, idx)
end

function _split_helper(::Type{<:AbstractArray}, buf_v, ::Val{N}, raw, idx) where {N}
    map(b -> _split_helper(eltype(b), b, Val(N - 1), raw, idx), buf_v)
end

function _split_helper(
        ::Type{<:AbstractArray}, buf_v::BlockedArray, ::Val{N}, raw, idx) where {N}
    BlockedArray(map(b -> _split_helper(eltype(b), b, Val(N - 1), raw, idx), buf_v),
        blocksizes(buf_v, 1))
end

function _split_helper(::Type{<:AbstractArray}, buf_v::Tuple, ::Val{N}, raw, idx) where {N}
    ntuple(i -> _split_helper(eltype(buf_v[i]), buf_v[i], Val(N - 1), raw, idx),
        Val(length(buf_v)))
end

function _split_helper(::Type{<:AbstractArray}, buf_v, ::Val{0}, raw, idx)
    _split_helper((), buf_v, (), raw, idx)
end

function _split_helper(_, buf_v, _, raw, idx)
    res = reshape(raw[idx[]:(idx[] + length(buf_v) - 1)], size(buf_v))
    idx[] += length(buf_v)
    return res
end

function _split_helper(_, buf_v::BlockedArray, _, raw, idx)
    res = BlockedArray(
        reshape(raw[idx[]:(idx[] + length(buf_v) - 1)], size(buf_v)), blocksizes(buf_v, 1))
    idx[] += length(buf_v)
    return res
end

function split_into_buffers(raw::AbstractArray, buf, recurse = Val(1))
    idx = Ref(1)
    ntuple(i -> _split_helper(buf[i], recurse, raw, idx), Val(length(buf)))
end

function _update_tuple_helper(buf_v::T, raw, idx) where {T}
    _update_tuple_helper(eltype(T), buf_v, raw, idx)
end

function _update_tuple_helper(::Type{<:AbstractArray}, buf_v, raw, idx)
    ntuple(i -> _update_tuple_helper(buf_v[i], raw, idx), length(buf_v))
end

function _update_tuple_helper(::Any, buf_v, raw, idx)
    copyto!(buf_v, view(raw, idx[]:(idx[] + length(buf_v) - 1)))
    idx[] += length(buf_v)
    return nothing
end

function update_tuple_of_buffers(raw::AbstractArray, buf)
    idx = Ref(1)
    ntuple(i -> _update_tuple_helper(buf[i], raw, idx), Val(length(buf)))
end

SciMLStructures.isscimlstructure(::MTKParameters) = true

SciMLStructures.ismutablescimlstructure(::MTKParameters) = true

function SciMLStructures.canonicalize(::SciMLStructures.Tunable, p::MTKParameters)
    arr = p.tunable
    repack = let p = p
        function (new_val)
            return SciMLStructures.replace(SciMLStructures.Tunable(), p, new_val)
        end
    end
    return arr, repack, true
end

function SciMLStructures.replace(::SciMLStructures.Tunable, p::MTKParameters, newvals)
    @set! p.tunable = newvals
    return p
end

function SciMLStructures.replace!(::SciMLStructures.Tunable, p::MTKParameters, newvals)
    copyto!(p.tunable, newvals)
    return nothing
end

function SciMLStructures.canonicalize(::SciMLStructures.Initials, p::MTKParameters)
    arr = p.initials
    repack = let p = p
        function (new_val)
            return SciMLStructures.replace(SciMLStructures.Initials(), p, new_val)
        end
    end
    return arr, repack, true
end

function SciMLStructures.replace(::SciMLStructures.Initials, p::MTKParameters, newvals)
    @set! p.initials = newvals
    return p
end

function SciMLStructures.replace!(::SciMLStructures.Initials, p::MTKParameters, newvals)
    copyto!(p.initials, newvals)
    return nothing
end

for (Portion, field, recurse) in [(SciMLStructures.Discrete, :discrete, 1)
     (SciMLStructures.Constants, :constant, 1)
     (Nonnumeric, :nonnumeric, 1)
     (SciMLStructures.Caches, :caches, 1)]
    @eval function SciMLStructures.canonicalize(::$Portion, p::MTKParameters)
        as_vector = buffer_to_arraypartition(p.$field)
        repack = let p = p
            function (new_val)
                return SciMLStructures.replace(($Portion)(), p, new_val)
            end
        end
        return as_vector, repack, true
    end

    @eval function SciMLStructures.replace(::$Portion, p::MTKParameters, newvals)
        @set! p.$field = split_into_buffers(newvals, p.$field, Val($recurse))
        p
    end

    @eval function SciMLStructures.replace!(::$Portion, p::MTKParameters, newvals)
        update_tuple_of_buffers(newvals, p.$field)
        nothing
    end
end

function Base.copy(p::MTKParameters)
    tunable = copy(p.tunable)
    initials = copy(p.initials)
    discrete = Tuple(eltype(buf) <: Real ? copy(buf) : copy.(buf) for buf in p.discrete)
    constant = Tuple(eltype(buf) <: Real ? copy(buf) : copy.(buf) for buf in p.constant)
    nonnumeric = isempty(p.nonnumeric) ? p.nonnumeric : copy.(p.nonnumeric)
    caches = isempty(p.caches) ? p.caches : copy.(p.caches)
    return MTKParameters(
        tunable,
        initials,
        discrete,
        constant,
        nonnumeric,
        caches
    )
end

function ArrayInterface.ismutable(::Type{MTKParameters{
        T, I, D, C, N, H}}) where {T, I, D, C, N, H}
    ArrayInterface.ismutable(T) || ArrayInterface.ismutable(I) ||
        any(ArrayInterface.ismutable, fieldtypes(D)) ||
        any(ArrayInterface.ismutable, fieldtypes(C)) ||
        any(ArrayInterface.ismutable, fieldtypes(N))
end

function SymbolicIndexingInterface.parameter_values(p::MTKParameters, pind::ParameterIndex)
    _ducktyped_parameter_values(p, pind)
end
function _ducktyped_parameter_values(p, pind::ParameterIndex)
    @unpack portion, idx = pind
    if portion isa SciMLStructures.Tunable
        return idx isa Int ? p.tunable[idx] : view(p.tunable, idx)
    end
    if portion isa SciMLStructures.Initials
        return idx isa Int ? p.initials[idx] : view(p.initials, idx)
    end
    i, j, k... = idx
    if portion isa SciMLStructures.Discrete
        return isempty(k) ? p.discrete[i][j] : p.discrete[i][j][k...]
    elseif portion isa SciMLStructures.Constants
        return isempty(k) ? p.constant[i][j] : p.constant[i][j][k...]
    elseif portion === NONNUMERIC_PORTION
        return isempty(k) ? p.nonnumeric[i][j] : p.nonnumeric[i][j][k...]
    else
        error("Unhandled portion $portion")
    end
end

function SymbolicIndexingInterface.set_parameter!(
        p::MTKParameters, val, pidx::ParameterIndex)
    @unpack portion, idx, validate_size = pidx
    if portion isa SciMLStructures.Tunable
        if validate_size && size(val) !== size(idx)
            throw(InvalidParameterSizeException(size(idx), size(val)))
        end
        p.tunable[idx] = val
    elseif portion isa SciMLStructures.Initials
        if validate_size && size(val) !== size(idx)
            throw(InvalidParameterSizeException(size(idx), size(val)))
        end
        p.initials[idx] = val
    else
        i, j, k... = idx
        if portion isa SciMLStructures.Discrete
            if isempty(k)
                if validate_size && size(val) !== size(p.discrete[i][j])
                    throw(InvalidParameterSizeException(
                        size(p.discrete[i][j]), size(val)))
                end
                p.discrete[i][j] = val
            else
                p.discrete[i][j][k...] = val
            end
        elseif portion isa SciMLStructures.Constants
            if isempty(k)
                if validate_size && size(val) !== size(p.constant[i][j])
                    throw(InvalidParameterSizeException(size(p.constant[i][j]), size(val)))
                end
                p.constant[i][j] = val
            else
                p.constant[i][j][k...] = val
            end
        elseif portion === NONNUMERIC_PORTION
            if isempty(k)
                p.nonnumeric[i][j] = val
            else
                p.nonnumeric[i][j][k...] = val
            end
        else
            error("Unhandled portion $portion")
        end
    end
    return nothing
end

function narrow_buffer_type_and_fallback_undefs(
        oldbuf::AbstractVector, newbuf::AbstractVector)
    type = Union{}
    for i in eachindex(newbuf)
        isassigned(newbuf, i) || continue
        type = promote_type(type, typeof(newbuf[i]))
    end
    if type == Union{}
        type = eltype(oldbuf)
    end
    newerbuf = similar(newbuf, type)
    for i in eachindex(newbuf)
        if isassigned(newbuf, i)
            newerbuf[i] = newbuf[i]
        else
            newerbuf[i] = oldbuf[i]
        end
    end
    return newerbuf
end

function validate_parameter_type(ic::IndexCache, p, idx::ParameterIndex, val)
    p = unwrap(p)
    if p isa Symbol
        p = get(ic.symbol_to_variable, p, nothing)
        p === nothing && return validate_parameter_type(ic, idx, val)
    end
    stype = symtype(p)
    sz = if stype <: AbstractArray
        Symbolics.shape(p) == Symbolics.Unknown() ? Symbolics.Unknown() : size(p)
    elseif stype <: Number
        size(p)
    else
        Symbolics.Unknown()
    end
    validate_parameter_type(ic, stype, sz, p, idx, val)
end

function validate_parameter_type(ic::IndexCache, idx::ParameterIndex, val)
    stype = get_buffer_template(ic, idx).type
    if (idx.portion == SciMLStructures.Tunable() ||
        idx.portion == SciMLStructures.Initials()) && !(idx.idx isa Int)
        stype = AbstractArray{<:stype}
    end
    validate_parameter_type(
        ic, stype, Symbolics.Unknown(), nothing, idx, val)
end

function validate_parameter_type(ic::IndexCache, stype, sz, sym, index, val)
    (; portion) = index
    if stype <: FnType
        stype = fntype_to_function_type(stype)
    end
    # Nonnumeric parameters have to match the type
    if portion === NONNUMERIC_PORTION
        val isa stype && return nothing
        throw(ParameterTypeException(
            :validate_parameter_type, sym === nothing ? index : sym, stype, val))
    end
    # Array parameters need array values...
    if stype <: AbstractArray && !isa(val, AbstractArray)
        throw(ParameterTypeException(
            :validate_parameter_type, sym === nothing ? index : sym, stype, val))
    end
    # ... and must match sizes
    if stype <: AbstractArray && sz != Symbolics.Unknown() && size(val) != sz
        throw(InvalidParameterSizeException(sym, val))
    end
    # Early exit
    val isa stype && return nothing
    if stype <: AbstractArray
        # Arrays need handling when eltype is `Real` (accept any real array)
        etype = eltype(stype)
        if etype <: Real
            etype = Real
        end
        # This is for duals and other complicated number types
        etype = SciMLBase.parameterless_type(etype)
        eltype(val) <: etype || throw(ParameterTypeException(
            :validate_parameter_type, sym === nothing ? index : sym, AbstractArray{etype}, val))
    else
        # Real check
        if stype <: Real
            stype = Real
        end
        stype = SciMLBase.parameterless_type(stype)
        val isa stype ||
            throw(ParameterTypeException(
                :validate_parameter_type, sym === nothing ? index : sym, stype, val))
    end
end

function indp_to_system(indp)
    while hasmethod(symbolic_container, Tuple{typeof(indp)})
        indp = symbolic_container(indp)
    end
    return indp
end

function SymbolicIndexingInterface.remake_buffer(indp, oldbuf::MTKParameters, idxs, vals)
    _remake_buffer(indp, oldbuf, idxs, vals)
end

# For type-inference when using `SII.setp_oop`
@generated function _remake_buffer(
        indp, oldbuf::MTKParameters{T, I, D, C, N, H},
        idxs::Union{Tuple{Vararg{ParameterIndex}}, AbstractArray{<:ParameterIndex{P}}},
        vals::Union{AbstractArray, Tuple}; validate = true) where {T, I, D, C, N, H, P}

    # fallback to non-generated method if values aren't type-stable
    if vals <: AbstractArray && !isconcretetype(eltype(vals))
        return quote
            $__remake_buffer(indp, oldbuf, collect(idxs), vals; validate)
        end
    end

    # given an index in idxs/vals and the current `eltype` of the buffer,
    # return the promoted eltype of the buffer
    function promote_valtype(i, valT)
        # tuples have distinct types, arrays have a common eltype
        valT′ = vals <: AbstractArray ? eltype(vals) : fieldtype(vals, i)
        # if the buffer is a scalarized buffer but the variable is an array
        # e.g. an array tunable, take the eltype
        if valT′ <: AbstractArray && !(valT <: AbstractArray)
            valT′ = eltype(valT′)
        end
        return promote_type(valT, valT′)
    end

    # types of the idxs
    idxtypes = if idxs <: AbstractArray
        # if both are arrays, there is only one possible type to check
        if vals <: AbstractArray
            (eltype(idxs),)
        else
            # if `vals` is a tuple, we repeat `eltype(idxs)` to check against
            # every possible type of the buffer
            ntuple(Returns(eltype(idxs)), Val(fieldcount(vals)))
        end
    else
        # `idxs` is a tuple, so we check against all buffers
        fieldtypes(idxs)
    end
    # promote types
    tunablesT = eltype(T)
    for (i, idxT) in enumerate(idxtypes)
        idxT <: ParameterIndex{SciMLStructures.Tunable} || continue
        tunablesT = promote_valtype(i, tunablesT)
    end
    initialsT = eltype(I)
    for (i, idxT) in enumerate(idxtypes)
        idxT <: ParameterIndex{SciMLStructures.Initials} || continue
        initialsT = promote_valtype(i, initialsT)
    end
    discretesT = ntuple(Val(fieldcount(D))) do i
        bufT = eltype(fieldtype(D, i))
        for (j, idxT) in enumerate(idxtypes)
            idxT <: ParameterIndex{SciMLStructures.Discrete, i} || continue
            bufT = promote_valtype(i, bufT)
        end
        bufT
    end
    constantsT = ntuple(Val(fieldcount(C))) do i
        bufT = eltype(fieldtype(C, i))
        for (j, idxT) in enumerate(idxtypes)
            idxT <: ParameterIndex{SciMLStructures.Constants, i} || continue
            bufT = promote_valtype(i, bufT)
        end
        bufT
    end
    nonnumericT = ntuple(Val(fieldcount(N))) do i
        bufT = eltype(fieldtype(N, i))
        for (j, idxT) in enumerate(idxtypes)
            idxT <: ParameterIndex{Nonnumeric, i} || continue
            bufT = promote_valtype(i, bufT)
        end
        bufT
    end

    expr = quote
        tunables = $similar(oldbuf.tunable, $tunablesT)
        copyto!(tunables, oldbuf.tunable)
        initials = $similar(oldbuf.initials, $initialsT)
        copyto!(initials, oldbuf.initials)
        discretes = $(Expr(:tuple,
            (:($similar(oldbuf.discrete[$i], $(discretesT[i])))
            for i in 1:length(discretesT))...))
        $((:($copyto!(discretes[$i], oldbuf.discrete[$i]))
        for i in 1:length(discretesT))...)
        constants = $(Expr(:tuple,
            (:($similar(oldbuf.constant[$i], $(constantsT[i])))
            for i in 1:length(constantsT))...))
        $((:($copyto!(constants[$i], oldbuf.constant[$i]))
        for i in 1:length(constantsT))...)
        nonnumerics = $(Expr(:tuple,
            (:($similar(oldbuf.nonnumeric[$i], $(nonnumericT[i])))
            for i in 1:length(nonnumericT))...))
        $((:($copyto!(nonnumerics[$i], oldbuf.nonnumeric[$i]))
        for i in 1:length(nonnumericT))...)
        caches = copy.(oldbuf.caches)
        newbuf = MTKParameters(
            tunables, initials, discretes, constants, nonnumerics, caches)
    end
    if idxs <: AbstractArray
        push!(expr.args, :(for (idx, val) in zip(idxs, vals)
            $setindex!(newbuf, val, idx)
        end))
    else
        for i in 1:fieldcount(idxs)
            push!(expr.args, :($setindex!(newbuf, vals[$i], idxs[$i])))
        end
    end
    if !ArrayInterface.ismutable(oldbuf)
        push!(expr.args, :(tunables = $similar_type($T, $tunablesT)(tunables)))
        push!(expr.args, :(initials = $similar_type($I, $initialsT)(initials)))
        push!(expr.args,
            :(discretes = $(Expr(:tuple,
                (:($similar_type($(fieldtype(D, i)), $(discretesT[i]))(discretes[$i]))
                for i in 1:length(discretesT))...))))
        push!(expr.args,
            :(constants = $(Expr(:tuple,
                (:($similar_type($(fieldtype(C, i)), $(constantsT[i]))(constants[$i]))
                for i in 1:length(constantsT))...))))
        push!(expr.args,
            :(nonnumerics = $(Expr(:tuple,
                (:($similar_type($(fieldtype(C, i)), $(nonnumericT[i]))(nonnumerics[$i]))
                for i in 1:length(nonnumericT))...))))
        push!(expr.args,
            :(newbuf = MTKParameters(
                tunables, initials, discretes, constants, nonnumerics, caches)))
    end
    push!(expr.args, :(return newbuf))

    return expr
end

function _remake_buffer(indp, oldbuf::MTKParameters, idxs, vals; validate = true)
    return __remake_buffer(indp, oldbuf, idxs, vals; validate)
end

function __remake_buffer(indp, oldbuf::MTKParameters, idxs, vals; validate = true)
    newbuf = @set oldbuf.tunable = similar(oldbuf.tunable, Any)
    @set! newbuf.initials = similar(oldbuf.initials, Any)
    @set! newbuf.discrete = Tuple(similar(buf, Any) for buf in newbuf.discrete)
    @set! newbuf.constant = Tuple(similar(buf, Any) for buf in newbuf.constant)
    @set! newbuf.nonnumeric = Tuple(similar(buf, Any) for buf in newbuf.nonnumeric)

    function handle_parameter(ic, sym, idx, val)
        if validate
            if sym === nothing
                validate_parameter_type(ic, idx, val)
            else
                validate_parameter_type(ic, sym, idx, val)
            end
        end
        # `ParameterIndex(idx)` turns off size validation since it relies on there
        # being an existing value
        set_parameter!(newbuf, val, ParameterIndex(idx))
    end

    handled_idxs = Set{ParameterIndex}()
    # If the parameter buffer is an `MTKParameters` object, `indp` must eventually drill
    # down to an `AbstractSystem` using `symbolic_container`. We leverage this to get
    # the index cache.
    ic = get_index_cache(indp_to_system(indp))
    for (idx, val) in zip(idxs, vals)
        sym = nothing
        if val === missing
            val = get_temporary_value(idx)
        end
        if symbolic_type(idx) == ScalarSymbolic()
            sym = idx
            idx = parameter_index(ic, sym)
            if idx === nothing
                @warn "Symbolic variable $sym is not a (non-dependent) parameter in the system"
                continue
            end
            idx in handled_idxs && continue
            handle_parameter(ic, sym, idx, val)
            push!(handled_idxs, idx)
        elseif symbolic_type(idx) == ArraySymbolic()
            sym = idx
            idx = parameter_index(ic, sym)
            if idx === nothing
                Symbolics.shape(sym) == Symbolics.Unknown() &&
                    throw(ParameterNotInSystem(sym))
                size(sym) == size(val) || throw(InvalidParameterSizeException(sym, val))

                for (i, vali) in zip(eachindex(sym), eachindex(val))
                    idx = parameter_index(ic, sym[i])
                    if idx === nothing
                        @warn "Symbolic variable $sym is not a (non-dependent) parameter in the system"
                        continue
                    end
                    # Intentionally don't check handled_idxs here because array variables always take priority
                    # See Issue#2804
                    handle_parameter(ic, sym[i], idx, val[vali])
                    push!(handled_idxs, idx)
                end
            else
                idx in handled_idxs && continue
                handle_parameter(ic, sym, idx, val)
                push!(handled_idxs, idx)
            end
        else # NotSymbolic
            if !(idx isa ParameterIndex)
                throw(ArgumentError("Expected index for parameter to be a symbolic variable or `ParameterIndex`, got $idx"))
            end
            handle_parameter(ic, nothing, idx, val)
        end
    end

    @set! newbuf.tunable = narrow_buffer_type_and_fallback_undefs(
        oldbuf.tunable, newbuf.tunable)
    if eltype(newbuf.tunable) <: Integer
        T = promote_type(eltype(newbuf.tunable), Float64)
        @set! newbuf.tunable = T.(newbuf.tunable)
    end
    @set! newbuf.initials = narrow_buffer_type_and_fallback_undefs(
        oldbuf.initials, newbuf.initials)
    if eltype(newbuf.initials) <: Integer
        T = promote_type(eltype(newbuf.initials), Float64)
        @set! newbuf.initials = T.(newbuf.initials)
    end
    @set! newbuf.discrete = narrow_buffer_type_and_fallback_undefs.(
        oldbuf.discrete, newbuf.discrete)
    @set! newbuf.constant = narrow_buffer_type_and_fallback_undefs.(
        oldbuf.constant, newbuf.constant)
    for (oldv, newv) in zip(oldbuf.nonnumeric, newbuf.nonnumeric)
        for i in eachindex(oldv)
            isassigned(newv, i) && continue
            newv[i] = oldv[i]
        end
    end
    @set! newbuf.nonnumeric = Tuple(
        typeof(oldv)(newv) for (oldv, newv) in zip(oldbuf.nonnumeric, newbuf.nonnumeric))
    if !ArrayInterface.ismutable(oldbuf)
        @set! newbuf.tunable = similar_type(oldbuf.tunable, eltype(newbuf.tunable))(newbuf.tunable)
        @set! newbuf.initials = similar_type(oldbuf.initials, eltype(newbuf.initials))(newbuf.initials)
        @set! newbuf.discrete = ntuple(Val(length(newbuf.discrete))) do i
            similar_type.(oldbuf.discrete[i], eltype(newbuf.discrete[i]))(newbuf.discrete[i])
        end
        @set! newbuf.constant = ntuple(Val(length(newbuf.constant))) do i
            similar_type.(oldbuf.constant[i], eltype(newbuf.constant[i]))(newbuf.constant[i])
        end
        @set! newbuf.nonnumeric = ntuple(Val(length(newbuf.nonnumeric))) do i
            similar_type.(oldbuf.nonnumeric[i], eltype(newbuf.nonnumeric[i]))(newbuf.nonnumeric[i])
        end
    end
    return newbuf
end

function as_any_buffer(p::MTKParameters)
    @set! p.tunable = similar(p.tunable, Any)
    @set! p.initials = similar(p.initials, Any)
    @set! p.discrete = Tuple(similar(buf, Any) for buf in p.discrete)
    @set! p.constant = Tuple(similar(buf, Any) for buf in p.constant)
    @set! p.nonnumeric = Tuple(similar(buf, Any) for buf in p.nonnumeric)
    @set! p.caches = Tuple(similar(buf, Any) for buf in p.caches)
    return p
end

struct NestedGetIndex{T}
    x::T
end

function Base.getindex(ngi::NestedGetIndex, idx::Tuple)
    i, j, k... = idx
    return ngi.x[i][j][k...]
end
function Base.getindex(ngi::NestedGetIndex, idx::NTuple{2})
    i, j = idx
    return ngi.x[i][j]
end

# Required for DiffEqArray constructor to work during interpolation
Base.size(::NestedGetIndex) = ()

function SymbolicIndexingInterface.with_updated_parameter_timeseries_values(
        ::AbstractSystem, ps::MTKParameters, args::Pair{A, B}...) where {
        A, B <: NestedGetIndex}
    for (i, ngi) in args
        for (j, val) in enumerate(ngi.x)
            copyto!(view(ps.discrete[j], Block(i)), val)
        end
    end
    return ps
end

function SciMLBase.create_parameter_timeseries_collection(
        sys::AbstractSystem, ps::MTKParameters, tspan)
    ic = get_index_cache(sys) # this exists because the parameters are `MTKParameters`
    isempty(ps.discrete) && return nothing
    num_discretes = only(blocksize(ps.discrete[1]))
    buffers = []
    partition_type = typeof(SciMLBase.get_saveable_values(sys, ps, 1))
    for i in 1:num_discretes
        ts = eltype(tspan)[]
        us = partition_type[]
        push!(buffers, DiffEqArray(us, ts, (1, 1)))
    end

    return ParameterTimeseriesCollection(Tuple(buffers), copy(ps))
end

@inline __get_blocks(tsidx::Int) = ()
@inline function __get_blocks(tsidx::Int, buffer::BlockedArray, buffers...)
    (buffer[Block(tsidx)], __get_blocks(tsidx, buffers...)...)
end
@inline function __get_blocks(tsidx::Int, buffer::BlockedArray{<:AbstractArray}, buffers...)
    (copy.(buffer[Block(tsidx)]), __get_blocks(tsidx, buffers...)...)
end

function SciMLBase.get_saveable_values(
        sys::AbstractSystem, ps::MTKParameters, timeseries_idx)
    return NestedGetIndex(__get_blocks(timeseries_idx, ps.discrete...))
end

function save_callback_discretes!(integ::SciMLBase.DEIntegrator, callback)
    ic = get_index_cache(indp_to_system(integ))
    ic === nothing && return
    clockidxs = get(ic.callback_to_clocks, callback, nothing)
    clockidxs === nothing && return

    for idx in clockidxs
        SciMLBase.save_discretes!(integ, idx)
    end
end

function DiffEqBase.anyeltypedual(
        p::MTKParameters, ::Type{Val{counter}} = Val{0}) where {counter}
    DiffEqBase.anyeltypedual(p.tunable)
end
function DiffEqBase.anyeltypedual(p::Type{<:MTKParameters{T}},
        ::Type{Val{counter}} = Val{0}) where {counter} where {T}
    DiffEqBase.__anyeltypedual(T)
end

# for compiling callbacks
# getindex indexes the vectors, setindex! linearly indexes values
# it's inconsistent, but we need it to be this way
@generated function Base.getindex(
        ps::MTKParameters{T, I, D, C, N, H}, idx::Int) where {T, I, D, C, N, H}
    paths = []
    if !(T <: SizedVector{0})
        push!(paths, :(ps.tunable))
    end
    if !(I <: SizedVector{0})
        push!(paths, :(ps.initials))
    end
    for i in 1:fieldcount(D)
        push!(paths, :(ps.discrete[$i]))
    end
    for i in 1:fieldcount(C)
        push!(paths, :(ps.constant[$i]))
    end
    for i in 1:fieldcount(N)
        push!(paths, :(ps.nonnumeric[$i]))
    end
    for i in 1:fieldcount(H)
        push!(paths, :(ps.caches[$i]))
    end
    expr = Expr(:if, :(idx == 1), :(return $(paths[1])))
    curexpr = expr
    for i in 2:length(paths)
        push!(curexpr.args, Expr(:elseif, :(idx == $i), :(return $(paths[i]))))
        curexpr = curexpr.args[end]
    end
    return Expr(:block, expr, :(throw(BoundsError(ps, idx))))
end

@generated function Base.length(ps::MTKParameters{
        T, I, D, C, N, H}) where {T, I, D, C, N, H}
    len = 0
    if !(T <: SizedVector{0})
        len += 1
    end
    if !(I <: SizedVector{0})
        len += 1
    end
    len += fieldcount(D) + fieldcount(C) + fieldcount(N) + fieldcount(H)
    return len
end

Base.size(ps::MTKParameters) = (length(ps),)

Base.IndexStyle(::Type{T}) where {T <: MTKParameters} = IndexLinear()

Base.getindex(p::MTKParameters, pind::ParameterIndex) = parameter_values(p, pind)

Base.setindex!(p::MTKParameters, val, pind::ParameterIndex) = set_parameter!(p, val, pind)

function Base.iterate(buf::MTKParameters, state = 1)
    total_len = length(buf)
    if state <= total_len
        return (buf[state], state + 1)
    else
        return nothing
    end
end

function Base.:(==)(a::MTKParameters, b::MTKParameters)
    return a.tunable == b.tunable && a.initials == b.initials && a.discrete == b.discrete &&
           a.constant == b.constant && a.nonnumeric == b.nonnumeric &&
           all(Iterators.map(a.caches, b.caches) do acache, bcache
               eltype(acache) == eltype(bcache) && length(acache) == length(bcache)
           end)
end

const MISSING_PARAMETERS_MESSAGE = """
                                Some parameters are missing from the variable map.
                                Please provide a value or default for the following variables:
                                """

struct MissingParametersError <: Exception
    vars::Any
end

function Base.showerror(io::IO, e::MissingParametersError)
    println(io, MISSING_PARAMETERS_MESSAGE)
    println(io, join(e.vars, ", "))
end

function InvalidParameterSizeException(param, val)
    DimensionMismatch("InvalidParameterSizeException: For parameter $(param) expected value of size $(size(param)). Received value $(val) of size $(size(val)).")
end

function InvalidParameterSizeException(param::Tuple, val::Tuple)
    DimensionMismatch("InvalidParameterSizeException: Expected value of size $(param). Received value of size $(val).")
end

function ParameterTypeException(func, param, expected, val)
    TypeError(func, "Parameter $param", expected, val)
end

struct ParameterNotInSystem <: Exception
    p::Any
end

function Base.showerror(io::IO, e::ParameterNotInSystem)
    println(io, "Symbolic variable $(e.p) is not a parameter in the system")
end
