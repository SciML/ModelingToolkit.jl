abstract type SymbolHash end

getsymbolhash(sym) = hasmetadata(sym, SymbolHash) ? getmetadata(sym, SymbolHash) : hash(unwrap(sym))

struct BufferTemplate
    type::DataType
    length::Int
end

struct ParameterIndex{P}
    portion::P
    idx::Int
end

struct IndexCache
    unknown_idx::Dict{UInt, Int}
    discrete_idx::Dict{UInt, Int}
    param_idx::Dict{UInt, Int}
    constant_idx::Dict{UInt, Int}
    dependent_idx::Dict{UInt, Int}
    discrete_buffer_sizes::Vector{BufferTemplate}
    param_buffer_sizes::Vector{BufferTemplate}
    constant_buffer_sizes::Vector{BufferTemplate}
    dependent_buffer_sizes::Vector{BufferTemplate}
end

function IndexCache(sys::AbstractSystem)
    unks = solved_unknowns(sys)
    unk_idxs = Dict{UInt, Int}()
    for (i, sym) in enumerate(unks)
        h = hash(unwrap(sym))
        unk_idxs[h] = i
        setmetadata(sym, SymbolHash, h)
    end

    disc_buffers = Dict{DataType, Set{BasicSymbolic}}()
    tunable_buffers = Dict{DataType, Set{BasicSymbolic}}()
    constant_buffers = Dict{DataType, Set{BasicSymbolic}}()
    dependent_buffers = Dict{DataType, Set{BasicSymbolic}}()

    function insert_by_type!(buffers::Dict{DataType, Set{BasicSymbolic}}, sym)
        sym = unwrap(sym)
        ctype = concrete_symtype(sym)
        buf = get!(buffers, ctype, Set{BasicSymbolic}())
        push!(buf, sym)
    end

    affs = vcat(affects(continuous_events(sys)), affects(discrete_events(sys)))
    for affect in affs
        if affect isa Equation
            is_parameter(sys, affect.lhs) || continue
            insert_by_type!(disc_buffers, affect.lhs)
        else
            discs = discretes(affect)
            for disc in discs
                is_parameter(sys, disc) || error("Expected discrete variable $disc in callback to be a parameter")
                insert_by_type!(disc_buffers, disc)
            end
        end
    end

    all_ps = Set(unwrap.(parameters(sys)))
    for (sym, value) in defaults(sys)
        sym = unwrap(sym)
        if sym in all_ps && symbolic_type(unwrap(value)) !== NotSymbolic()
            insert_by_type!(dependent_buffers, sym)
        end
    end

    for p in parameters(sys)
        p = unwrap(p)
        ctype = concrete_symtype(p)
        haskey(disc_buffers, ctype) && p in disc_buffers[ctype] && continue
        haskey(dependent_buffers, ctype) && p in dependent_buffers[ctype] && continue

        insert_by_type!(
            if is_discrete_domain(p)
                disc_buffers
            elseif istunable(p, true)
                tunable_buffers
            else
                constant_buffers
            end,
            p
        )
    end

    function get_buffer_sizes_and_idxs(buffers::Dict{DataType, Set{BasicSymbolic}})
        idxs = Dict{UInt, Int}()
        buffer_sizes = BufferTemplate[]
        idx = 1
        for (T, buf) in buffers
            for p in buf
                h = hash(p)
                setmetadata(p, SymbolHash, h)
                idxs[h] = idx
                idx += 1
            end
            push!(buffer_sizes, BufferTemplate(T, length(buf)))
        end
        return idxs, buffer_sizes
    end

    disc_idxs, discrete_buffer_sizes = get_buffer_sizes_and_idxs(disc_buffers)
    param_idxs, param_buffer_sizes = get_buffer_sizes_and_idxs(tunable_buffers)
    const_idxs, const_buffer_sizes = get_buffer_sizes_and_idxs(constant_buffers)
    dependent_idxs, dependent_buffer_sizes = get_buffer_sizes_and_idxs(dependent_buffers)

    return IndexCache(
        unk_idxs,
        disc_idxs,
        param_idxs,
        const_idxs,
        dependent_idxs,
        discrete_buffer_sizes,
        param_buffer_sizes,
        const_buffer_sizes,
        dependent_buffer_sizes,
    )
end

function reorder_parameters(sys::AbstractSystem, ps; kwargs...)
    if has_index_cache(sys) && get_index_cache(sys) !== nothing
        reorder_parameters(get_index_cache(sys), ps; kwargs...)
    elseif ps isa Tuple
        ps
    else
        (ps,)
    end
end

function reorder_parameters(ic::IndexCache, ps; drop_missing = false)
    param_buf = ArrayPartition((fill(variable(:DEF), temp.length) for temp in ic.param_buffer_sizes)...)
    disc_buf = ArrayPartition((fill(variable(:DEF), temp.length) for temp in ic.discrete_buffer_sizes)...)
    const_buf = ArrayPartition((fill(variable(:DEF), temp.length) for temp in ic.constant_buffer_sizes)...)
    dep_buf = ArrayPartition((fill(variable(:DEF), temp.length) for temp in ic.dependent_buffer_sizes)...)

    for p in ps
        h = getsymbolhash(p)
        if haskey(ic.discrete_idx, h)
            disc_buf[ic.discrete_idx[h]] = unwrap(p)
        elseif haskey(ic.param_idx, h)
            param_buf[ic.param_idx[h]] = unwrap(p)
        elseif haskey(ic.constant_idx, h)
            const_buf[ic.constant_idx[h]] = unwrap(p)
        elseif haskey(ic.dependent_idx, h)
            dep_buf[ic.dependent_idx[h]] = unwrap(p)
        else
            error("Invalid parameter $p")
        end
    end

    result = broadcast.(unwrap, (param_buf.x..., disc_buf.x..., const_buf.x..., dep_buf.x...))
    if drop_missing
        result = map(result) do buf
            filter(buf) do sym
                return !isequal(sym, unwrap(variable(:DEF)))
            end
        end
    end
    return result
end

concrete_symtype(x::BasicSymbolic) = concrete_symtype(symtype(x))
concrete_symtype(::Type{Real}) = Float64
concrete_symtype(::Type{Integer}) = Int
concrete_symtype(::Type{Vector{T}}) where {T} = Vector{concrete_symtype(T)}
concrete_symtype(::Type{T}) where {T} = T
