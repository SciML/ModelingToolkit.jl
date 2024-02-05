abstract type SymbolHash end

struct IndexCache
    unknown_idx::Dict{UInt, Int}
    discrete_idx::Dict{UInt, Int}
    param_idx::Dict{UInt, Int}
    discrete_buffer_type_and_size::Vector{Tuple{DataType, Int}}
    param_buffer_type_and_size::Vector{Tuple{DataType, Int}}
end

function IndexCache(sys::AbstractSystem)
    unks = solved_unknowns(sys)
    unk_idxs = Dict{UInt, Int}()
    for (i, sym) in enumerate(unks)
        h = hash(unwrap(sym))
        unk_idxs[h] = i
        setmetadata(sym, SymbolHash, h)
    end

    # split parameters, also by type
    discrete_params = Dict{DataType, Any}()
    tunable_params = Dict{DataType, Any}()

    for p in parameters(sys)
        T = symtype(p)
        buf = get!(is_discrete_domain(p) ? discrete_params : tunable_params, T, [])
        push!(buf, unwrap(p))
    end
    
    disc_idxs = Dict{UInt, Int}()
    discrete_buffer_type_and_size = Tuple{DataType, Int}[]
    didx = 1

    for (T, ps) in discrete_params
        push!(discrete_buffer_type_and_size, (T, length(ps)))
        for p in ps
            h = hash(p)
            disc_idxs[h] = didx
            didx += 1
            setmetadata(p, SymbolHash, h)
        end
    end

    param_idxs = Dict{UInt, Int}()
    param_buffer_type_and_size = Tuple{DataType, Int}[]
    pidx = 1

    for (T, ps) in tunable_params
        push!(param_buffer_type_and_size, (T, length(ps)))
        for p in ps
            h = hash(p)
            param_idxs[h] = pidx
            pidx += 1
            setmetadata(p, SymbolHash, h)
        end
    end

    return IndexCache(unk_idxs, disc_idxs, param_idxs, discrete_buffer_type_and_size, param_buffer_type_and_size)
end

function reorder_parameters(ic::IndexCache, ps)
    param_bufs = ArrayPartition((Vector{BasicSymbolic{T}}(undef, sz) for (T, sz) in ic.param_buffer_type_and_size)...)
    disc_bufs = ArrayPartition((Vector{BasicSymbolic{T}}(undef, sz) for (T, sz) in ic.discrete_buffer_type_and_size)...)

    for p in ps
        h = hasmetadata(p, SymbolHash) ? getmetadata(p, SymbolHash) : hash(unwrap(p))
        if haskey(ic.discrete_idx, h)
            disc_bufs[ic.discrete_idx[h]] = unwrap(p)
        elseif haskey(ic.param_idx, h)
            param_bufs[ic.param_idx[h]] = unwrap(p)
        else
            error("Invalid parameter $p")
        end
    end

    return (param_bufs.x..., disc_bufs.x...)
end
