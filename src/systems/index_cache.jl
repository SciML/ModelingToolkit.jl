abstract type SymbolHash end

function getsymbolhash(sym)
    sym = unwrap(sym)
    hasmetadata(sym, SymbolHash) ? getmetadata(sym, SymbolHash) : hash(sym)
end

struct BufferTemplate
    type::DataType
    length::Int
end

const DEPENDENT_PORTION = :dependent
const NONNUMERIC_PORTION = :nonnumeric

struct ParameterIndex{P, I}
    portion::P
    idx::I
end

const IndexMap = Dict{UInt, Tuple{Int, Int}}

struct IndexCache
    unknown_idx::Dict{UInt, Int}
    discrete_idx::IndexMap
    param_idx::IndexMap
    constant_idx::IndexMap
    dependent_idx::IndexMap
    nonnumeric_idx::IndexMap
    discrete_buffer_sizes::Vector{BufferTemplate}
    param_buffer_sizes::Vector{BufferTemplate}
    constant_buffer_sizes::Vector{BufferTemplate}
    dependent_buffer_sizes::Vector{BufferTemplate}
    nonnumeric_buffer_sizes::Vector{BufferTemplate}
end

function IndexCache(sys::AbstractSystem)
    unks = solved_unknowns(sys)
    unk_idxs = Dict{UInt, Int}()
    for (i, sym) in enumerate(unks)
        h = getsymbolhash(sym)
        unk_idxs[h] = i
    end

    disc_buffers = Dict{DataType, Set{BasicSymbolic}}()
    tunable_buffers = Dict{DataType, Set{BasicSymbolic}}()
    constant_buffers = Dict{DataType, Set{BasicSymbolic}}()
    dependent_buffers = Dict{DataType, Set{BasicSymbolic}}()
    nonnumeric_buffers = Dict{DataType, Set{BasicSymbolic}}()

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
                is_parameter(sys, disc) ||
                    error("Expected discrete variable $disc in callback to be a parameter")
                insert_by_type!(disc_buffers, disc)
            end
        end
    end
    if has_discrete_subsystems(sys) && get_discrete_subsystems(sys) !== nothing
        _, inputs, continuous_id, _ = get_discrete_subsystems(sys)
        for par in inputs[continuous_id]
            is_parameter(sys, par) || error("Discrete subsytem input is not a parameter")
            insert_by_type!(disc_buffers, par)
        end
    end

    # all_ps = Set(unwrap.(parameters(sys)))
    # for (sym, value) in defaults(sys)
    #     sym = unwrap(sym)
    #     if sym in all_ps && symbolic_type(unwrap(value)) !== NotSymbolic()
    #         insert_by_type!(dependent_buffers, sym)
    #     end
    # end

    for p in parameters(sys)
        p = unwrap(p)
        ctype = concrete_symtype(p)
        haskey(disc_buffers, ctype) && p in disc_buffers[ctype] && continue
        haskey(dependent_buffers, ctype) && p in dependent_buffers[ctype] && continue

        insert_by_type!(
            if ctype <: Real || ctype <: Vector{<:Real}
                if is_discrete_domain(p)
                    disc_buffers
                elseif istunable(p, true) && size(p) !== Symbolics.Unknown()
                    tunable_buffers
                else
                    constant_buffers
                end
            else
                nonnumeric_buffers
            end,
            p
        )
    end

    function get_buffer_sizes_and_idxs(
            buffers::Dict{DataType, Set{BasicSymbolic}}, track_linear_index = true)
        idxs = IndexMap()
        buffer_sizes = BufferTemplate[]
        for (i, (T, buf)) in enumerate(buffers)
            for (j, p) in enumerate(buf)
                h = getsymbolhash(p)
                idxs[h] = (i, j)
                h = getsymbolhash(default_toterm(p))
                idxs[h] = (i, j)
            end
            push!(buffer_sizes, BufferTemplate(T, length(buf)))
        end
        return idxs, buffer_sizes
    end

    disc_idxs, discrete_buffer_sizes = get_buffer_sizes_and_idxs(disc_buffers)
    param_idxs, param_buffer_sizes = get_buffer_sizes_and_idxs(tunable_buffers)
    const_idxs, const_buffer_sizes = get_buffer_sizes_and_idxs(constant_buffers)
    dependent_idxs, dependent_buffer_sizes = get_buffer_sizes_and_idxs(dependent_buffers)
    nonnumeric_idxs, nonnumeric_buffer_sizes = get_buffer_sizes_and_idxs(nonnumeric_buffers)

    return IndexCache(
        unk_idxs,
        disc_idxs,
        param_idxs,
        const_idxs,
        dependent_idxs,
        nonnumeric_idxs,
        discrete_buffer_sizes,
        param_buffer_sizes,
        const_buffer_sizes,
        dependent_buffer_sizes,
        nonnumeric_buffer_sizes
    )
end

function ParameterIndex(ic::IndexCache, p)
    p = unwrap(p)
    if istree(p) && operation(p) === getindex
        sub_idx = Base.tail(arguments(p))
        p = arguments(p)[begin]
    else
        sub_idx = ()
    end
    h = getsymbolhash(p)
    return if haskey(ic.param_idx, h)
        ParameterIndex(SciMLStructures.Tunable(), (ic.param_idx[h]..., sub_idx...))
    elseif haskey(ic.discrete_idx, h)
        ParameterIndex(SciMLStructures.Discrete(), (ic.discrete_idx[h]..., sub_idx...))
    elseif haskey(ic.constant_idx, h)
        ParameterIndex(SciMLStructures.Constants(), (ic.constant_idx[h]..., sub_idx...))
    elseif haskey(ic.dependent_idx, h)
        ParameterIndex(DEPENDENT_PORTION, (ic.dependent_idx[h]..., sub_idx...))
    elseif haskey(ic.nonnumeric_idx, h)
        ParameterIndex(NONNUMERIC_PORTION, (ic.nonnumeric_idx[h]..., sub_idx...))
    else
        nothing
    end
end

function discrete_linear_index(ic::IndexCache, idx::ParameterIndex)
    idx.portion isa SciMLStructures.Discrete || error("Discrete variable index expected")
    ind = sum(temp.length for temp in ic.param_buffer_sizes; init = 0)
    ind += sum(
        temp.length for temp in Iterators.take(ic.discrete_buffer_sizes, idx.idx[1] - 1);
        init = 0)
    ind += idx.idx[2]
    return ind
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
    param_buf = Tuple(BasicSymbolic[unwrap(variable(:DEF)) for _ in 1:(temp.length)]
    for temp in ic.param_buffer_sizes)
    disc_buf = Tuple(BasicSymbolic[unwrap(variable(:DEF)) for _ in 1:(temp.length)]
    for temp in ic.discrete_buffer_sizes)
    const_buf = Tuple(BasicSymbolic[unwrap(variable(:DEF)) for _ in 1:(temp.length)]
    for temp in ic.constant_buffer_sizes)
    dep_buf = Tuple(BasicSymbolic[unwrap(variable(:DEF)) for _ in 1:(temp.length)]
    for temp in ic.dependent_buffer_sizes)
    nonnumeric_buf = Tuple(BasicSymbolic[unwrap(variable(:DEF)) for _ in 1:(temp.length)]
    for temp in ic.nonnumeric_buffer_sizes)

    for p in ps
        h = getsymbolhash(p)
        if haskey(ic.discrete_idx, h)
            i, j = ic.discrete_idx[h]
            disc_buf[i][j] = unwrap(p)
        elseif haskey(ic.param_idx, h)
            i, j = ic.param_idx[h]
            param_buf[i][j] = unwrap(p)
        elseif haskey(ic.constant_idx, h)
            i, j = ic.constant_idx[h]
            const_buf[i][j] = unwrap(p)
        elseif haskey(ic.dependent_idx, h)
            i, j = ic.dependent_idx[h]
            dep_buf[i][j] = unwrap(p)
        elseif haskey(ic.nonnumeric_idx, h)
            i, j = ic.nonnumeric_idx[h]
            nonnumeric_buf[i][j] = unwrap(p)
        else
            error("Invalid parameter $p")
        end
    end

    result = broadcast.(
        unwrap, (param_buf..., disc_buf..., const_buf..., nonnumeric_buf..., dep_buf...))
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
