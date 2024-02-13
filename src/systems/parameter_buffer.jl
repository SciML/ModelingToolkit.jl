struct MTKParameters{T, D, C, E, F}
    tunable::T
    discrete::D
    constant::C
    dependent::E
    dependent_update::F
end

function MTKParameters(sys::AbstractSystem, p; tofloat = false, use_union = false)
    ic = if has_index_cache(sys) && get_index_cache(sys) !== nothing
        get_index_cache(sys)
    else
        error("Cannot create MTKParameters if system does not have index_cache")
    end
    all_ps = Set(unwrap.(parameters(sys)))
    union!(all_ps, default_toterm.(unwrap.(parameters(sys))))
    if p isa Vector && !(eltype(p) <: Pair) && !isempty(p)
        ps = parameters(sys)
        length(p) == length(ps) || error("Invalid parameters")
        p = ps .=> p
    end
    defs = Dict(default_toterm(unwrap(k)) => v for (k, v) in defaults(sys) if unwrap(k) in all_ps)
    if p isa SciMLBase.NullParameters
        p = defs
    else
        extra_params = Dict(unwrap(k) => v for (k, v) in p if !in(unwrap(k), all_ps))
        p = merge(defs, Dict(default_toterm(unwrap(k)) => v for (k, v) in p if unwrap(k) in all_ps))
        p = Dict(k => fixpoint_sub(v, extra_params) for (k, v) in p if !haskey(extra_params, unwrap(k)))
    end

    tunable_buffer = ArrayPartition((Vector{temp.type}(undef, temp.length) for temp in ic.param_buffer_sizes)...)
    disc_buffer = ArrayPartition((Vector{temp.type}(undef, temp.length) for temp in ic.discrete_buffer_sizes)...)
    const_buffer = ArrayPartition((Vector{temp.type}(undef, temp.length) for temp in ic.constant_buffer_sizes)...)
    dep_buffer = ArrayPartition((Vector{temp.type}(undef, temp.length) for temp in ic.dependent_buffer_sizes)...)
    dependencies = Dict{Num, Num}()
    function set_value(sym, val)
        h = getsymbolhash(sym)
        if haskey(ic.param_idx, h)
            tunable_buffer[ic.param_idx[h]] = val
        elseif haskey(ic.discrete_idx, h)
            disc_buffer[ic.discrete_idx[h]] = val
        elseif haskey(ic.constant_idx, h)
            const_buffer[ic.constant_idx[h]] = val
        elseif haskey(ic.dependent_idx, h)
            dep_buffer[ic.dependent_idx[h]] = val
            dependencies[wrap(sym)] = wrap(p[sym])
        elseif !isequal(default_toterm(sym), sym)
            set_value(default_toterm(sym), val)
        else
            error("Symbol $sym does not have an index")
        end
    end

    for (sym, val) in p
        sym = unwrap(sym)
        ctype = concrete_symtype(sym)
        val = convert(ctype, fixpoint_sub(val, p))
        if size(sym) == ()
            set_value(sym, val)
        else
            if length(sym) != length(val)
                error("Size of $sym does not match size of initial value $val")
            end
            for (i, j) in zip(eachindex(sym), eachindex(val))
                set_value(sym[i], val[j])
            end
        end
    end

    dep_exprs = ArrayPartition((wrap.(v) for v in dep_buffer.x)...)
    for (sym, val) in dependencies
        h = getsymbolhash(sym)
        idx = ic.dependent_idx[h]
        dep_exprs[idx] = wrap(fixpoint_sub(val, dependencies))
    end
    p = reorder_parameters(ic, parameters(sys))[begin:end-length(dep_buffer.x)]
    update_function = if isempty(dep_exprs.x)
        (_...) -> ()
    else
        RuntimeGeneratedFunctions.@RuntimeGeneratedFunction(build_function(dep_exprs, p...)[2])
    end
    # everything is an ArrayPartition so it's easy to figure out how many
    # distinct vectors we have for each portion as `ArrayPartition.x`
    if isempty(tunable_buffer.x)
        tunable_buffer = ArrayPartition(Float64[])
    end
    if isempty(disc_buffer.x)
        disc_buffer = ArrayPartition(Float64[])
    end
    if isempty(const_buffer.x)
        const_buffer = ArrayPartition(Float64[])
    end
    if isempty(dep_buffer.x)
        dep_buffer = ArrayPartition(Float64[])
    end
    if use_union
        tunable_buffer = ArrayPartition(restrict_array_to_union(tunable_buffer))
        disc_buffer = ArrayPartition(restrict_array_to_union(disc_buffer))
        const_buffer = ArrayPartition(restrict_array_to_union(const_buffer))
        dep_buffer = ArrayPartition(restrict_array_to_union(dep_buffer))
    elseif tofloat
        tunable_buffer = ArrayPartition(Float64.(tunable_buffer))
        disc_buffer = ArrayPartition(Float64.(disc_buffer))
        const_buffer = ArrayPartition(Float64.(const_buffer))
        dep_buffer = ArrayPartition(Float64.(dep_buffer))
    end
    return MTKParameters{typeof(tunable_buffer), typeof(disc_buffer), typeof(const_buffer),
                         typeof(dep_buffer), typeof(update_function)}(tunable_buffer,
                         disc_buffer, const_buffer, dep_buffer, update_function)
end

SciMLStructures.isscimlstructure(::MTKParameters) = true

SciMLStructures.ismutablescimlstructure(::MTKParameters) = true

for (Portion, field) in [
    (SciMLStructures.Tunable, :tunable)
    (SciMLStructures.Discrete, :discrete)
    (SciMLStructures.Constants, :constant)
]
    @eval function SciMLStructures.canonicalize(::$Portion, p::MTKParameters)
        function repack(values)
            p.$field .= values
            p.dependent_update(p.dependent, p.tunable.x..., p.discrete.x..., p.constant.x...)
        end
        return p.$field, repack, true
    end

    @eval function SciMLStructures.replace(::$Portion, p::MTKParameters, newvals)
        @set! p.$field = newvals
        p.dependent_update(p.dependent, p.tunable.x..., p.discrete.x..., p.constant.x...)
    end

    @eval function SciMLStructures.replace!(::$Portion, p::MTKParameters, newvals)
        p.$field .= newvals
        p.dependent_update(p.dependent, p.tunable.x..., p.discrete.x..., p.constant.x...)
        nothing
    end
end

function SymbolicIndexingInterface.parameter_values(p::MTKParameters, i::ParameterIndex)
    @unpack portion, idx = i
    if portion isa SciMLStructures.Tunable
        return p.tunable[idx]
    elseif portion isa SciMLStructures.Discrete
        return p.discrete[idx]
    elseif portion isa SciMLStructures.Constants
        return p.constant[idx]
    elseif portion === nothing
        return p.dependent[idx]
    else
        error("Unhandled portion $portion")
    end
end

function SymbolicIndexingInterface.set_parameter!(p::MTKParameters, val, idx::ParameterIndex)
    @unpack portion, idx = idx
    if portion isa SciMLStructures.Tunable
        p.tunable[idx] = val
    elseif portion isa SciMLStructures.Discrete
        p.discrete[idx] = val
    elseif portion isa SciMLStructures.Constants
        p.constant[idx] = val
    elseif portion === nothing
        error("Cannot set value of parameter: ")
    else
        error("Unhandled portion $portion")
    end
    p.dependent_update(p.dependent, p.tunable.x..., p.discrete.x..., p.constant.x...)
end

# for compiling callbacks
# getindex indexes the vectors, setindex! linearly indexes values
# it's inconsistent, but we need it to be this way
function Base.getindex(buf::MTKParameters, i)
    if !isempty(buf.tunable)
        i <= length(buf.tunable.x) && return buf.tunable.x[i]
        i -= length(buf.tunable.x)
    end
    if !isempty(buf.discrete)
        i <= length(buf.discrete.x) && return buf.discrete.x[i]
        i -= length(buf.discrete.x)
    end
    if !isempty(buf.constant)
        i <= length(buf.constant.x) && return buf.constant.x[i]
        i -= length(buf.constant.x)
    end
    isempty(buf.dependent) || return buf.dependent.x[i]
    throw(BoundsError(buf, i))
end
function Base.setindex!(buf::MTKParameters, val, i)
    if i <= length(buf.tunable)
        buf.tunable[i] = val
    elseif i <= length(buf.tunable) + length(buf.discrete)
        buf.discrete[i - length(buf.tunable)] = val
    else
        buf.constant[i - length(buf.tunable) - length(buf.discrete)] = val
    end
    buf.dependent_update(buf.dependent, buf.tunable.x..., buf.discrete.x..., buf.constant.x...)
end

function Base.iterate(buf::MTKParameters, state = 1)
    total_len = 0
    isempty(buf.tunable) || (total_len += length(buf.tunable.x))
    isempty(buf.discrete) || (total_len += length(buf.discrete.x))
    isempty(buf.constant) || (total_len += length(buf.constant.x))
    isempty(buf.dependent) || (total_len += length(buf.dependent.x))
    if state <= total_len
        return (buf[state], state + 1)
    else
        return nothing
    end
end

function Base.:(==)(a::MTKParameters, b::MTKParameters)
    return a.tunable == b.tunable && a.discrete == b.discrete &&
        a.constant == b.constant && a.dependent == b.dependent
end

# to support linearize/linearization_function
function jacobian_wrt_vars(pf::F, p::MTKParameters, input_idxs, chunk::C) where {F, C}
    T = eltype(p.tunable)
    tag = ForwardDiff.Tag(pf, T)
    dualtype = ForwardDiff.Dual{typeof(tag), T, ForwardDiff.chunksize(chunk)}
    tunable, _, _ = SciMLStructures.canonicalize(SciMLStructures.Tunable(), p)
    p_big = SciMLStructures.replace(SciMLStructures.Tunable(), p, dualtype.(tunable))
    p_closure = let pf = pf,
        input_idxs = input_idxs,
        p_big = p_big

        function (p_small_inner)
            tunable, repack, _ = SciMLStructures.canonicalize(SciMLStructures.Tunable(), p_big)
            tunable[input_idxs] .= p_small_inner
            p_big = repack(tunable)
            pf(p_big)
        end
    end
    tunable, _, _ = SciMLStructures.canonicalize(SciMLStructures.Tunable(), p)
    p_small = tunable[input_idxs]
    cfg = ForwardDiff.JacobianConfig(p_closure, p_small, chunk, tag)
    ForwardDiff.jacobian(p_closure, p_small, cfg, Val(false))
end

function as_duals(p::MTKParameters, dualtype)
    tunable = dualtype.(p.tunable)
    discrete = dualtype.(p.discrete)
    return MTKParameters{typeof(tunable), typeof(discrete)}(tunable, discrete)
end
