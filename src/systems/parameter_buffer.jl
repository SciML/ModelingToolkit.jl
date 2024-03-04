struct MTKParameters{T, D, C, E, N, F, G}
    tunable::T
    discrete::D
    constant::C
    dependent::E
    nonnumeric::N
    dependent_update_iip::F
    dependent_update_oop::G
end

function MTKParameters(sys::AbstractSystem, p; tofloat = false, use_union = false)
    ic = if has_index_cache(sys) && get_index_cache(sys) !== nothing
        get_index_cache(sys)
    else
        error("Cannot create MTKParameters if system does not have index_cache")
    end
    all_ps = Set(unwrap.(full_parameters(sys)))
    union!(all_ps, default_toterm.(unwrap.(full_parameters(sys))))
    if p isa Vector && !(eltype(p) <: Pair) && !isempty(p)
        ps = full_parameters(sys)
        length(p) == length(ps) || error("Invalid parameters")
        p = ps .=> p
    end
    defs = Dict(default_toterm(unwrap(k)) => v
    for (k, v) in defaults(sys)
    if unwrap(k) in all_ps || default_toterm(unwrap(k)) in all_ps)
    if p isa SciMLBase.NullParameters
        p = defs
    else
        extra_params = Dict(unwrap(k) => v
        for (k, v) in p if !in(unwrap(k), all_ps) && !in(default_toterm(unwrap(k)), all_ps))
        p = merge(defs,
            Dict(default_toterm(unwrap(k)) => v
            for (k, v) in p if unwrap(k) in all_ps || default_toterm(unwrap(k)) in all_ps))
        p = Dict(k => fixpoint_sub(v, extra_params)
        for (k, v) in p if !haskey(extra_params, unwrap(k)))
    end

    for (sym, _) in p
        if istree(sym) && operation(sym) === getindex &&
           first(arguments(sym)) in all_ps
            error("Scalarized parameter values ($sym) are not supported. Instead of `[p[1] => 1.0, p[2] => 2.0]` use `[p => [1.0, 2.0]]`")
        end
    end

    tunable_buffer = Tuple(Vector{temp.type}(undef, temp.length)
    for temp in ic.tunable_buffer_sizes)
    disc_buffer = Tuple(Vector{temp.type}(undef, temp.length)
    for temp in ic.discrete_buffer_sizes)
    const_buffer = Tuple(Vector{temp.type}(undef, temp.length)
    for temp in ic.constant_buffer_sizes)
    dep_buffer = Tuple(Vector{temp.type}(undef, temp.length)
    for temp in ic.dependent_buffer_sizes)
    nonnumeric_buffer = Tuple(Vector{temp.type}(undef, temp.length)
    for temp in ic.nonnumeric_buffer_sizes)
    function set_value(sym, val)
        done = true
        if haskey(ic.tunable_idx, sym)
            i, j = ic.tunable_idx[sym]
            tunable_buffer[i][j] = val
        elseif haskey(ic.discrete_idx, sym)
            i, j = ic.discrete_idx[sym]
            disc_buffer[i][j] = val
        elseif haskey(ic.constant_idx, sym)
            i, j = ic.constant_idx[sym]
            const_buffer[i][j] = val
        elseif haskey(ic.dependent_idx, sym)
            i, j = ic.dependent_idx[sym]
            dep_buffer[i][j] = val
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
        ctype = concrete_symtype(sym)
        val = convert(ctype, fixpoint_sub(val, p))
        done = set_value(sym, val)
        if !done && Symbolics.isarraysymbolic(sym)
            done = all(set_value.(collect(sym), val))
        end
        if !done
            error("Symbol $sym does not have an index")
        end
    end

    if has_parameter_dependencies(sys) &&
       (pdeps = get_parameter_dependencies(sys)) !== nothing
        pdeps = Dict(k => fixpoint_sub(v, pdeps) for (k, v) in pdeps)
        dep_exprs = ArrayPartition((wrap.(v) for v in dep_buffer)...)
        for (sym, val) in pdeps
            i, j = ic.dependent_idx[sym]
            dep_exprs.x[i][j] = wrap(val)
        end
        p = reorder_parameters(ic, full_parameters(sys))
        oop, iip = build_function(dep_exprs, p...)
        update_function_iip, update_function_oop = RuntimeGeneratedFunctions.@RuntimeGeneratedFunction(iip),
        RuntimeGeneratedFunctions.@RuntimeGeneratedFunction(oop)
    else
        update_function_iip = update_function_oop = nothing
    end

    mtkps = MTKParameters{
        typeof(tunable_buffer), typeof(disc_buffer), typeof(const_buffer),
        typeof(dep_buffer), typeof(nonnumeric_buffer), typeof(update_function_iip),
        typeof(update_function_oop)}(tunable_buffer, disc_buffer, const_buffer, dep_buffer,
        nonnumeric_buffer, update_function_iip, update_function_oop)
    if mtkps.dependent_update_iip !== nothing
        mtkps.dependent_update_iip(ArrayPartition(mtkps.dependent), mtkps...)
    end
    return mtkps
end

function buffer_to_arraypartition(buf)
    return ArrayPartition(Tuple(eltype(v) <: AbstractArray ? buffer_to_arraypartition(v) :
                                v for v in buf))
end

function split_into_buffers(raw::AbstractArray, buf; recurse = true)
    idx = 1
    function _helper(buf_v; recurse = true)
        if eltype(buf_v) <: AbstractArray && recurse
            return _helper.(buf_v; recurse = false)
        else
            res = reshape(raw[idx:(idx + length(buf_v) - 1)], size(buf_v))
            idx += length(buf_v)
            return res
        end
    end
    return Tuple(_helper(buf_v; recurse) for buf_v in buf)
end

function update_tuple_of_buffers(raw::AbstractArray, buf)
    idx = 1
    function _helper(buf_v)
        if eltype(buf_v) <: AbstractArray
            _helper.(buf_v)
        else
            copyto!(buf_v, view(raw, idx:(idx + length(buf_v) - 1)))
            idx += length(buf_v)
        end
    end
    _helper.(buf)
end

SciMLStructures.isscimlstructure(::MTKParameters) = true

SciMLStructures.ismutablescimlstructure(::MTKParameters) = true

for (Portion, field) in [(SciMLStructures.Tunable, :tunable)
                         (SciMLStructures.Discrete, :discrete)
                         (SciMLStructures.Constants, :constant)]
    @eval function SciMLStructures.canonicalize(::$Portion, p::MTKParameters)
        as_vector = buffer_to_arraypartition(p.$field)
        repack = let as_vector = as_vector, p = p
            function (new_val)
                if new_val !== as_vector
                    update_tuple_of_buffers(new_val, p.$field)
                end
                if p.dependent_update_iip !== nothing
                    p.dependent_update_iip(ArrayPartition(p.dependent), p...)
                end
                p
            end
        end
        return as_vector, repack, true
    end

    @eval function SciMLStructures.replace(::$Portion, p::MTKParameters, newvals)
        @set! p.$field = split_into_buffers(newvals, p.$field)
        if p.dependent_update_oop !== nothing
            raw = p.dependent_update_oop(p...)
            @set! p.dependent = split_into_buffers(raw, p.dependent; recurse = false)
        end
        p
    end

    @eval function SciMLStructures.replace!(::$Portion, p::MTKParameters, newvals)
        src = split_into_buffers(newvals, p.$field)
        for i in 1:length(p.$field)
            (p.$field)[i] .= src[i]
        end
        if p.dependent_update_iip !== nothing
            p.dependent_update_iip(ArrayPartition(p.dependent), p...)
        end
        nothing
    end
end

function Base.copy(p::MTKParameters)
    tunable = Tuple(eltype(buf) <: Real ? copy(buf) : copy.(buf) for buf in p.tunable)
    discrete = Tuple(eltype(buf) <: Real ? copy(buf) : copy.(buf) for buf in p.discrete)
    constant = Tuple(eltype(buf) <: Real ? copy(buf) : copy.(buf) for buf in p.constant)
    dependent = Tuple(eltype(buf) <: Real ? copy(buf) : copy.(buf) for buf in p.dependent)
    nonnumeric = copy.(p.nonnumeric)
    return MTKParameters(
        tunable,
        discrete,
        constant,
        dependent,
        nonnumeric,
        p.dependent_update_iip,
        p.dependent_update_oop
    )
end

function SymbolicIndexingInterface.parameter_values(p::MTKParameters, pind::ParameterIndex)
    @unpack portion, idx = pind
    i, j, k... = idx
    if portion isa SciMLStructures.Tunable
        return p.tunable[i][j][k...]
    elseif portion isa SciMLStructures.Discrete
        return p.discrete[i][j][k...]
    elseif portion isa SciMLStructures.Constants
        return p.constant[i][j][k...]
    elseif portion === DEPENDENT_PORTION
        return p.dependent[i][j][k...]
    elseif portion === NONNUMERIC_PORTION
        return isempty(k) ? p.nonnumeric[i][j] : p.nonnumeric[i][j][k...]
    else
        error("Unhandled portion $portion")
    end
end

function SymbolicIndexingInterface.set_parameter!(
        p::MTKParameters, val, idx::ParameterIndex)
    @unpack portion, idx = idx
    i, j, k... = idx
    if portion isa SciMLStructures.Tunable
        if isempty(k)
            p.tunable[i][j] = val
        else
            p.tunable[i][j][k...] = val
        end
    elseif portion isa SciMLStructures.Discrete
        if isempty(k)
            p.discrete[i][j] = val
        else
            p.discrete[i][j][k...] = val
        end
    elseif portion isa SciMLStructures.Constants
        if isempty(k)
            p.constant[i][j] = val
        else
            p.constant[i][j][k...] = val
        end
    elseif portion === DEPENDENT_PORTION
        error("Cannot set value of dependent parameter")
    elseif portion === NONNUMERIC_PORTION
        if isempty(k)
            p.nonnumeric[i][j] = val
        else
            p.nonnumeric[i][j][k...] = val
        end
    else
        error("Unhandled portion $portion")
    end
    if p.dependent_update_iip !== nothing
        p.dependent_update_iip(ArrayPartition(p.dependent), p...)
    end
end

function _set_parameter_unchecked!(
        p::MTKParameters, val, idx::ParameterIndex; update_dependent = true)
    @unpack portion, idx = idx
    i, j, k... = idx
    if portion isa SciMLStructures.Tunable
        if isempty(k)
            p.tunable[i][j] = val
        else
            p.tunable[i][j][k...] = val
        end
    elseif portion isa SciMLStructures.Discrete
        if isempty(k)
            p.discrete[i][j] = val
        else
            p.discrete[i][j][k...] = val
        end
    elseif portion isa SciMLStructures.Constants
        if isempty(k)
            p.constant[i][j] = val
        else
            p.constant[i][j][k...] = val
        end
    elseif portion === DEPENDENT_PORTION
        if isempty(k)
            p.dependent[i][j] = val
        else
            p.dependent[i][j][k...] = val
        end
        update_dependent = false
    elseif portion === NONNUMERIC_PORTION
        if isempty(k)
            p.nonnumeric[i][j] = val
        else
            p.nonnumeric[i][j][k...] = val
        end
    else
        error("Unhandled portion $portion")
    end
    update_dependent && p.dependent_update_iip !== nothing &&
        p.dependent_update_iip(ArrayPartition(p.dependent), p...)
end

_subarrays(v::AbstractVector) = isempty(v) ? () : (v,)
_subarrays(v::ArrayPartition) = v.x
_subarrays(v::Tuple) = v
_num_subarrays(v::AbstractVector) = 1
_num_subarrays(v::ArrayPartition) = length(v.x)
_num_subarrays(v::Tuple) = length(v)
# for compiling callbacks
# getindex indexes the vectors, setindex! linearly indexes values
# it's inconsistent, but we need it to be this way
function Base.getindex(buf::MTKParameters, i)
    if !isempty(buf.tunable)
        i <= _num_subarrays(buf.tunable) && return _subarrays(buf.tunable)[i]
        i -= _num_subarrays(buf.tunable)
    end
    if !isempty(buf.discrete)
        i <= _num_subarrays(buf.discrete) && return _subarrays(buf.discrete)[i]
        i -= _num_subarrays(buf.discrete)
    end
    if !isempty(buf.constant)
        i <= _num_subarrays(buf.constant) && return _subarrays(buf.constant)[i]
        i -= _num_subarrays(buf.constant)
    end
    if !isempty(buf.nonnumeric)
        i <= _num_subarrays(buf.nonnumeric) && return _subarrays(buf.nonnumeric)[i]
        i -= _num_subarrays(buf.nonnumeric)
    end
    if !isempty(buf.dependent)
        i <= _num_subarrays(buf.dependent) && return _subarrays(buf.dependent)[i]
        i -= _num_subarrays(buf.dependent)
    end
    throw(BoundsError(buf, i))
end
function Base.setindex!(p::MTKParameters, val, i)
    function _helper(buf)
        done = false
        for v in buf
            if i <= length(v)
                v[i] = val
                done = true
            else
                i -= length(v)
            end
        end
        done
    end
    _helper(p.tunable) || _helper(p.discrete) || _helper(p.constant) ||
        _helper(p.nonnumeric) || throw(BoundsError(p, i))
    if p.dependent_update_iip !== nothing
        p.dependent_update_iip(ArrayPartition(p.dependent), p...)
    end
end

function Base.iterate(buf::MTKParameters, state = 1)
    total_len = 0
    total_len += _num_subarrays(buf.tunable)
    total_len += _num_subarrays(buf.discrete)
    total_len += _num_subarrays(buf.constant)
    total_len += _num_subarrays(buf.nonnumeric)
    total_len += _num_subarrays(buf.dependent)
    if state <= total_len
        return (buf[state], state + 1)
    else
        return nothing
    end
end

function Base.:(==)(a::MTKParameters, b::MTKParameters)
    return a.tunable == b.tunable && a.discrete == b.discrete &&
           a.constant == b.constant && a.dependent == b.dependent &&
           a.nonnumeric == b.nonnumeric
end

# to support linearize/linearization_function
function jacobian_wrt_vars(pf::F, p::MTKParameters, input_idxs, chunk::C) where {F, C}
    tunable, _, _ = SciMLStructures.canonicalize(SciMLStructures.Tunable(), p)
    T = eltype(tunable)
    tag = ForwardDiff.Tag(pf, T)
    dualtype = ForwardDiff.Dual{typeof(tag), T, ForwardDiff.chunksize(chunk)}
    p_big = SciMLStructures.replace(SciMLStructures.Tunable(), p, dualtype.(tunable))
    p_closure = let pf = pf,
        input_idxs = input_idxs,
        p_big = p_big

        function (p_small_inner)
            for (i, val) in zip(input_idxs, p_small_inner)
                _set_parameter_unchecked!(p_big, val, i)
            end
            # tunable, repack, _ = SciMLStructures.canonicalize(SciMLStructures.Tunable(), p_big)
            # tunable[input_idxs] .= p_small_inner
            # p_big = repack(tunable)
            return if pf isa SciMLBase.ParamJacobianWrapper
                buffer = Array{dualtype}(undef, size(pf.u))
                pf(buffer, p_big)
                buffer
            else
                pf(p_big)
            end
        end
    end
    # tunable, _, _ = SciMLStructures.canonicalize(SciMLStructures.Tunable(), p)
    # p_small = tunable[input_idxs]
    p_small = parameter_values.((p,), input_idxs)
    cfg = ForwardDiff.JacobianConfig(p_closure, p_small, chunk, tag)
    ForwardDiff.jacobian(p_closure, p_small, cfg, Val(false))
end

function as_duals(p::MTKParameters, dualtype)
    tunable = dualtype.(p.tunable)
    discrete = dualtype.(p.discrete)
    return MTKParameters{typeof(tunable), typeof(discrete)}(tunable, discrete)
end
