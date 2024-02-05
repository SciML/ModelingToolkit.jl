struct MTKParameters{T, D}
    tunable::T
    discrete::D
end

function MTKParameters(sys::AbstractSystem, p; toterm = default_toterm)
    ic = if has_index_cache(sys)
        get_index_cache(sys)
    else
        IndexCache(sys)
    end
    tunable_buffer = if length(ic.param_buffer_type_and_size) == 0
        Float64[]
    elseif length(ic.param_buffer_type_and_size) == 1
        T, sz = only(ic.param_buffer_type_and_size)
        Vector{T == Real ? Float64 : T}(undef, sz)
    else
        ArrayPartition((Vector{T == Real ? Float64 : T}(undef, sz) for (T, sz) in ic.param_buffer_type_and_size)...)
    end

    disc_buffer = if length(ic.discrete_buffer_type_and_size) == 0
        Float64[]
    elseif length(ic.discrete_buffer_type_and_size) == 1
        T, sz = only(ic.discrete_buffer_type_and_size)
        Vector{T == Real ? Float64 : T}(undef, sz)
    else
        ArrayPartition((Vector{T == Real ? Float64 : T}(undef, sz) for (T, sz) in ic.discrete_buffer_type_and_size)...)
    end

    for (sym, value) in defaults(sys)
        sym = toterm(unwrap(sym))
        h = hasmetadata(sym, SymbolHash) ? getmetadata(sym, SymbolHash) : hash(sym)
        if haskey(ic.discrete_idx, h)
            disc_buffer[ic.discrete_idx[h]] = value
        elseif haskey(ic.param_idx, h)
            tunable_buffer[ic.param_idx[h]] = value
        end
    end

    if !isa(p, SciMLBase.NullParameters)
        for (sym, value) in p
            sym = toterm(unwrap(sym))
            h = hasmetadata(sym, SymbolHash) ? getmetadata(sym, SymbolHash) : hash(sym)
            if haskey(ic.discrete_idx, h)
                disc_buffer[ic.discrete_idx[h]] = value
            elseif haskey(ic.param_idx, h)
                tunable_buffer[ic.param_idx[h]] = value
            else
                error("Invalid parameter $sym")
            end
        end
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
        return p.$field, repack, !isa(p.$field, ArrayPartition)
    end

    @eval function SciMLStructures.replace(::$Portion, p::MTKParameters, newvals)
        new_field = similar(p.$field)
        new_field .= newvals
        @set p.$field = new_field
    end

    @eval function SciMLStructures.replace!(::$Portion, p::MTKParameters, newvals)
        p.$field .= newvals
        nothing
    end
end

function raw_vectors(buf::MTKParameters)
    tunable = if isempty(buf.tunable)
        ()
    elseif buf.tunable isa ArrayPartition
        buf.tunable.x
    else
        (buf.tunable,)
    end
    discrete = if isempty(buf.discrete)
        ()
    elseif buf.discrete isa ArrayPartition
        buf.discrete.x
    else
        (buf.discrete,)
    end
    return (tunable..., discrete...)
end
