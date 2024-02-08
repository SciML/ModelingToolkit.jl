struct MTKParameters{T, D}
    tunable::T
    discrete::D
end

function MTKParameters(sys::AbstractSystem, p; toterm = default_toterm, tofloat = false, use_union = false)
    ic = if has_index_cache(sys) && get_index_cache(sys) !== nothing
        get_index_cache(sys)
    else
        error("Cannot create MTKParameters if system does not have index_cache")
    end
    all_ps = Set(unwrap.(parameters(sys)))
    if p isa Vector && !(eltype(p) <: Pair)
        ps = parameters(sys)
        length(p) == length(ps) || error("Invalid parameters")
        p = ps .=> p
    end
    defs = Dict(unwrap(k) => v for (k, v) in defaults(sys) if unwrap(k) in all_ps)
    if p isa SciMLBase.NullParameters
        p = defs
    else
        p = merge(defs, Dict(unwrap(k) => v for (k, v) in p))
    end

    tunable_buffer = ArrayPartition((Vector{temp.type}(undef, temp.length) for temp in ic.param_buffer_sizes)...)
    disc_buffer = ArrayPartition((Vector{temp.type}(undef, temp.length) for temp in ic.discrete_buffer_sizes)...)
    function set_value(sym, val)
        h = getsymbolhash(sym)
        if haskey(ic.param_idx, h)
            tunable_buffer[ic.param_idx[h]] = val
        elseif haskey(ic.discrete_idx, h)
            disc_buffer[ic.discrete_idx[h]] = val
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

    # everything is an ArrayPartition so it's easy to figure out how many
    # distinct vectors we have for each portion as `ArrayPartition.x`
    if tunable_buffer isa ArrayPartition && isempty(tunable_buffer.x) || isempty(tunable_buffer)
        tunable_buffer = ArrayPartition(Float64[])
    end
    if disc_buffer isa ArrayPartition && isempty(disc_buffer.x) || isempty(disc_buffer)
        disc_buffer = ArrayPartition(Float64[])
    end
    if use_union
        tunable_buffer = ArrayPartition(restrict_array_to_union(tunable_buffer))
        disc_buffer = ArrayPartition(restrict_array_to_union(disc_buffer))
    elseif tofloat
        tunable_buffer = ArrayPartition(Float64.(tunable_buffer))
        disc_buffer = ArrayPartition(Float64.(disc_buffer))
    end
    return MTKParameters{typeof(tunable_buffer), typeof(disc_buffer)}(tunable_buffer, disc_buffer)
end

SciMLStructures.isscimlstructure(::MTKParameters) = true

SciMLStructures.ismutablescimlstructure(::MTKParameters) = true

for (Portion, field) in [
    (SciMLStructures.Tunable, :tunable)
    (SciMLStructures.Discrete, :discrete)
]
    @eval function SciMLStructures.canonicalize(::$Portion, p::MTKParameters)
        function repack(values)
            p.$field .= values
        end
        return p.$field, repack, true
    end

    @eval function SciMLStructures.replace(::$Portion, p::MTKParameters, newvals)
        @set p.$field = newvals
    end

    @eval function SciMLStructures.replace!(::$Portion, p::MTKParameters, newvals)
        p.$field .= newvals
        nothing
    end
end

function SymbolicIndexingInterface.parameter_values(p::MTKParameters, i::ParameterIndex)
    @unpack portion, idx = i
    if portion isa SciMLStructures.Tunable
        return p.tunable[idx]
    elseif portion isa SciMLStructures.Discrete
        return p.discrete[idx]
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
    else
        error("Unhandled portion $portion")
    end
end

# for compiling callbacks
# getindex indexes the vectors, setindex! linearly indexes values
# it's inconsistent, but we need it to be this way
function Base.getindex(buf::MTKParameters, i)
    if i <= length(buf.tunable.x)
        buf.tunable.x[i]
    else
        buf.discrete.x[i - length(buf.tunable.x)]
    end
end
function Base.setindex!(buf::MTKParameters, val, i)
    if i <= length(buf.tunable)
        buf.tunable[i] = val
    else
        buf.discrete[i - length(buf.tunable)] = val
    end
end

function Base.iterate(buf::MTKParameters, state = 1)
    tunable = if isempty(buf.tunable)
        ()
    elseif buf.tunable isa ArrayPartition
        buf.tunable.x
    end
    discrete = if isempty(buf.discrete)
        ()
    elseif buf.discrete isa ArrayPartition
        buf.discrete.x
    end
    if state <= length(tunable)
        return (tunable[state], state + 1)
    elseif state <= length(tunable) + length(discrete)
        return (discrete[state - length(tunable)], state + 1)
    else
        return nothing
    end
end

function Base.:(==)(a::MTKParameters, b::MTKParameters)
    return a.tunable == b.tunable && a.discrete == b.discrete
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
