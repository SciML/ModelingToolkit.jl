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

const ParamIndexMap = Dict{Union{Symbol, BasicSymbolic}, Tuple{Int, Int}}
const UnknownIndexMap = Dict{Union{Symbol, BasicSymbolic}, Union{Int, UnitRange{Int}}}

struct IndexCache
    unknown_idx::UnknownIndexMap
    discrete_idx::ParamIndexMap
    tunable_idx::ParamIndexMap
    constant_idx::ParamIndexMap
    dependent_idx::ParamIndexMap
    nonnumeric_idx::ParamIndexMap
    discrete_buffer_sizes::Vector{BufferTemplate}
    tunable_buffer_sizes::Vector{BufferTemplate}
    constant_buffer_sizes::Vector{BufferTemplate}
    dependent_buffer_sizes::Vector{BufferTemplate}
    nonnumeric_buffer_sizes::Vector{BufferTemplate}
end

function IndexCache(sys::AbstractSystem)
    unks = solved_unknowns(sys)
    unk_idxs = UnknownIndexMap()
    let idx = 1
        for sym in unks
            usym = unwrap(sym)
            sym_idx = if Symbolics.isarraysymbolic(sym)
                idx:(idx + length(sym) - 1)
            else
                idx
            end
            unk_idxs[usym] = sym_idx

            if hasname(sym)
                unk_idxs[getname(usym)] = sym_idx
            end
            idx += length(sym)
        end
    end

    disc_buffers = Dict{Any, Set{BasicSymbolic}}()
    tunable_buffers = Dict{Any, Set{BasicSymbolic}}()
    constant_buffers = Dict{Any, Set{BasicSymbolic}}()
    dependent_buffers = Dict{Any, Set{BasicSymbolic}}()
    nonnumeric_buffers = Dict{Any, Set{BasicSymbolic}}()

    function insert_by_type!(buffers::Dict{Any, Set{BasicSymbolic}}, sym)
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
            is_parameter(sys, par) || error("Discrete subsystem input is not a parameter")
            insert_by_type!(disc_buffers, par)
        end
    end

    if has_parameter_dependencies(sys) &&
       (pdeps = get_parameter_dependencies(sys)) !== nothing
        for (sym, value) in pdeps
            sym = unwrap(sym)
            insert_by_type!(dependent_buffers, sym)
        end
    end

    for p in parameters(sys)
        p = unwrap(p)
        ctype = concrete_symtype(p)
        haskey(disc_buffers, ctype) && p in disc_buffers[ctype] && continue
        haskey(dependent_buffers, ctype) && p in dependent_buffers[ctype] && continue
        insert_by_type!(
            if ctype <: Real || ctype <: AbstractArray{<:Real}
                if is_discrete_domain(p)
                    disc_buffers
                elseif istunable(p, true) && Symbolics.shape(p) !== Symbolics.Unknown()
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

    function get_buffer_sizes_and_idxs(buffers::Dict{Any, Set{BasicSymbolic}})
        idxs = ParamIndexMap()
        buffer_sizes = BufferTemplate[]
        for (i, (T, buf)) in enumerate(buffers)
            for (j, p) in enumerate(buf)
                idxs[p] = (i, j)
                idxs[default_toterm(p)] = (i, j)
                if hasname(p)
                    idxs[getname(p)] = (i, j)
                    idxs[getname(default_toterm(p))] = (i, j)
                end
            end
            push!(buffer_sizes, BufferTemplate(T, length(buf)))
        end
        return idxs, buffer_sizes
    end

    disc_idxs, discrete_buffer_sizes = get_buffer_sizes_and_idxs(disc_buffers)
    tunable_idxs, tunable_buffer_sizes = get_buffer_sizes_and_idxs(tunable_buffers)
    const_idxs, const_buffer_sizes = get_buffer_sizes_and_idxs(constant_buffers)
    dependent_idxs, dependent_buffer_sizes = get_buffer_sizes_and_idxs(dependent_buffers)
    nonnumeric_idxs, nonnumeric_buffer_sizes = get_buffer_sizes_and_idxs(nonnumeric_buffers)

    return IndexCache(
        unk_idxs,
        disc_idxs,
        tunable_idxs,
        const_idxs,
        dependent_idxs,
        nonnumeric_idxs,
        discrete_buffer_sizes,
        tunable_buffer_sizes,
        const_buffer_sizes,
        dependent_buffer_sizes,
        nonnumeric_buffer_sizes
    )
end

function SymbolicIndexingInterface.is_variable(ic::IndexCache, sym)
    return check_index_map(ic.unknown_idx, sym) !== nothing
end

function SymbolicIndexingInterface.variable_index(ic::IndexCache, sym)
    return check_index_map(ic.unknown_idx, sym)
end

function SymbolicIndexingInterface.is_parameter(ic::IndexCache, sym)
    return check_index_map(ic.tunable_idx, sym) !== nothing ||
           check_index_map(ic.discrete_idx, sym) !== nothing ||
           check_index_map(ic.constant_idx, sym) !== nothing ||
           check_index_map(ic.nonnumeric_idx, sym) !== nothing ||
           check_index_map(ic.dependent_idx, sym) !== nothing
end

function SymbolicIndexingInterface.parameter_index(ic::IndexCache, sym)
    return if (idx = check_index_map(ic.tunable_idx, sym)) !== nothing
        ParameterIndex(SciMLStructures.Tunable(), idx)
    elseif (idx = check_index_map(ic.discrete_idx, sym)) !== nothing
        ParameterIndex(SciMLStructures.Discrete(), idx)
    elseif (idx = check_index_map(ic.constant_idx, sym)) !== nothing
        ParameterIndex(SciMLStructures.Constants(), idx)
    elseif (idx = check_index_map(ic.nonnumeric_idx, sym)) !== nothing
        ParameterIndex(NONNUMERIC_PORTION, idx)
    elseif (idx = check_index_map(ic.dependent_idx, sym)) !== nothing
        ParameterIndex(DEPENDENT_PORTION, idx)
    else
        nothing
    end
end

function check_index_map(idxmap, sym)
    if (idx = get(idxmap, sym, nothing)) !== nothing
        return idx
    elseif hasname(sym) && (idx = get(idxmap, getname(sym), nothing)) !== nothing
        return idx
    end
    dsym = default_toterm(sym)
    isequal(sym, dsym) && return nothing
    if (idx = get(idxmap, dsym, nothing)) !== nothing
        idx
    elseif hasname(dsym) && (idx = get(idxmap, getname(dsym), nothing)) !== nothing
        idx
    else
        nothing
    end
end

function ParameterIndex(ic::IndexCache, p, sub_idx = ())
    p = unwrap(p)
    return if haskey(ic.tunable_idx, p)
        ParameterIndex(SciMLStructures.Tunable(), (ic.tunable_idx[p]..., sub_idx...))
    elseif haskey(ic.discrete_idx, p)
        ParameterIndex(SciMLStructures.Discrete(), (ic.discrete_idx[p]..., sub_idx...))
    elseif haskey(ic.constant_idx, p)
        ParameterIndex(SciMLStructures.Constants(), (ic.constant_idx[p]..., sub_idx...))
    elseif haskey(ic.dependent_idx, p)
        ParameterIndex(DEPENDENT_PORTION, (ic.dependent_idx[p]..., sub_idx...))
    elseif haskey(ic.nonnumeric_idx, p)
        ParameterIndex(NONNUMERIC_PORTION, (ic.nonnumeric_idx[p]..., sub_idx...))
    elseif istree(p) && operation(p) === getindex
        _p, sub_idx... = arguments(p)
        ParameterIndex(ic, _p, sub_idx)
    else
        nothing
    end
end

function discrete_linear_index(ic::IndexCache, idx::ParameterIndex)
    idx.portion isa SciMLStructures.Discrete || error("Discrete variable index expected")
    ind = sum(temp.length for temp in ic.tunable_buffer_sizes; init = 0)
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
    isempty(ps) && return ()
    param_buf = Tuple(BasicSymbolic[unwrap(variable(:DEF)) for _ in 1:(temp.length)]
    for temp in ic.tunable_buffer_sizes)
    disc_buf = Tuple(BasicSymbolic[unwrap(variable(:DEF)) for _ in 1:(temp.length)]
    for temp in ic.discrete_buffer_sizes)
    const_buf = Tuple(BasicSymbolic[unwrap(variable(:DEF)) for _ in 1:(temp.length)]
    for temp in ic.constant_buffer_sizes)
    dep_buf = Tuple(BasicSymbolic[unwrap(variable(:DEF)) for _ in 1:(temp.length)]
    for temp in ic.dependent_buffer_sizes)
    nonnumeric_buf = Tuple(BasicSymbolic[unwrap(variable(:DEF)) for _ in 1:(temp.length)]
    for temp in ic.nonnumeric_buffer_sizes)

    for p in ps
        if haskey(ic.discrete_idx, p)
            i, j = ic.discrete_idx[p]
            disc_buf[i][j] = unwrap(p)
        elseif haskey(ic.tunable_idx, p)
            i, j = ic.tunable_idx[p]
            param_buf[i][j] = unwrap(p)
        elseif haskey(ic.constant_idx, p)
            i, j = ic.constant_idx[p]
            const_buf[i][j] = unwrap(p)
        elseif haskey(ic.dependent_idx, p)
            i, j = ic.dependent_idx[p]
            dep_buf[i][j] = unwrap(p)
        elseif haskey(ic.nonnumeric_idx, p)
            i, j = ic.nonnumeric_idx[p]
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
    if all(isempty, result)
        return ()
    end
    return result
end

concrete_symtype(x::BasicSymbolic) = concrete_symtype(symtype(x))
concrete_symtype(::Type{Real}) = Float64
concrete_symtype(::Type{Integer}) = Int
concrete_symtype(::Type{A}) where {T, N, A <: Array{T, N}} = Array{concrete_symtype(T), N}
concrete_symtype(::Type{T}) where {T} = T
