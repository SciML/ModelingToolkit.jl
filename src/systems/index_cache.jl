struct BufferTemplate
    type::Union{DataType, UnionAll}
    length::Int
end

function BufferTemplate(s::Type{<:Symbolics.Struct}, length::Int)
    T = Symbolics.juliatype(s)
    BufferTemplate(T, length)
end

struct Dependent <: SciMLStructures.AbstractPortion end
struct Nonnumeric <: SciMLStructures.AbstractPortion end
const DEPENDENT_PORTION = Dependent()
const NONNUMERIC_PORTION = Nonnumeric()

struct ParameterIndex{P, I}
    portion::P
    idx::I
    validate_size::Bool
end

ParameterIndex(portion, idx) = ParameterIndex(portion, idx, false)

const ParamIndexMap = Dict{BasicSymbolic, Tuple{Int, Int}}
const UnknownIndexMap = Dict{
    BasicSymbolic, Union{Int, UnitRange{Int}, AbstractArray{Int}}}

struct IndexCache
    unknown_idx::UnknownIndexMap
    discrete_idx::Dict{BasicSymbolic, Tuple{Int, Int, Int}}
    tunable_idx::ParamIndexMap
    constant_idx::ParamIndexMap
    dependent_idx::ParamIndexMap
    nonnumeric_idx::ParamIndexMap
    observed_syms::Set{BasicSymbolic}
    discrete_buffer_sizes::Vector{Vector{BufferTemplate}}
    tunable_buffer_sizes::Vector{BufferTemplate}
    constant_buffer_sizes::Vector{BufferTemplate}
    dependent_buffer_sizes::Vector{BufferTemplate}
    nonnumeric_buffer_sizes::Vector{BufferTemplate}
    symbol_to_variable::Dict{Symbol, BasicSymbolic}
end

function IndexCache(sys::AbstractSystem)
    unks = solved_unknowns(sys)
    unk_idxs = UnknownIndexMap()
    symbol_to_variable = Dict{Symbol, BasicSymbolic}()

    let idx = 1
        for sym in unks
            usym = unwrap(sym)
            rsym = renamespace(sys, usym)
            sym_idx = if Symbolics.isarraysymbolic(sym)
                reshape(idx:(idx + length(sym) - 1), size(sym))
            else
                idx
            end
            unk_idxs[usym] = sym_idx
            unk_idxs[rsym] = sym_idx
            idx += length(sym)
        end
        for sym in unks
            usym = unwrap(sym)
            iscall(sym) && operation(sym) === getindex || continue
            arrsym = arguments(sym)[1]
            all(haskey(unk_idxs, arrsym[i]) for i in eachindex(arrsym)) || continue

            idxs = [unk_idxs[arrsym[i]] for i in eachindex(arrsym)]
            if idxs == idxs[begin]:idxs[end]
                idxs = reshape(idxs[begin]:idxs[end], size(idxs))
            end
            rsym = renamespace(sys, arrsym)
            unk_idxs[arrsym] = idxs
            unk_idxs[rsym] = idxs
        end
    end

    observed_syms = Set{Union{Symbol, BasicSymbolic}}()
    for eq in observed(sys)
        if symbolic_type(eq.lhs) != NotSymbolic()
            sym = eq.lhs
            ttsym = default_toterm(sym)
            rsym = renamespace(sys, sym)
            rttsym = renamespace(sys, ttsym)
            push!(observed_syms, sym)
            push!(observed_syms, ttsym)
            push!(observed_syms, rsym)
            push!(observed_syms, rttsym)
        end
    end

    disc_buffers = Dict{Int, Dict{Any, Set{BasicSymbolic}}}()
    disc_clocks = Dict{Union{Symbol, BasicSymbolic}, Int}()
    tunable_buffers = Dict{Any, Set{BasicSymbolic}}()
    constant_buffers = Dict{Any, Set{BasicSymbolic}}()
    dependent_buffers = Dict{Any, Set{BasicSymbolic}}()
    nonnumeric_buffers = Dict{Any, Set{BasicSymbolic}}()

    function insert_by_type!(buffers::Dict{Any, Set{BasicSymbolic}}, sym)
        sym = unwrap(sym)
        ctype = symtype(sym)
        buf = get!(buffers, ctype, Set{BasicSymbolic}())
        push!(buf, sym)
    end

    if has_discrete_subsystems(sys) && get_discrete_subsystems(sys) !== nothing
        syss, inputs, continuous_id, _ = get_discrete_subsystems(sys)

        for (i, (inps, disc_sys)) in enumerate(zip(inputs, syss))
            i == continuous_id && continue
            disc_buffers[i] = Dict{Any, Set{BasicSymbolic}}()

            for inp in inps
                inp = unwrap(inp)
                ttinp = default_toterm(inp)
                rinp = renamespace(sys, inp)
                rttinp = renamespace(sys, ttinp)
                is_parameter(sys, inp) ||
                    error("Discrete subsystem $i input $inp is not a parameter")

                disc_clocks[inp] = i
                disc_clocks[ttinp] = i
                disc_clocks[rinp] = i
                disc_clocks[rttinp] = i

                insert_by_type!(disc_buffers[i], inp)
            end

            for sym in unknowns(disc_sys)
                sym = unwrap(sym)
                ttsym = default_toterm(sym)
                rsym = renamespace(sys, sym)
                rttsym = renamespace(sys, ttsym)
                is_parameter(sys, sym) ||
                    error("Discrete subsystem $i unknown $sym is not a parameter")

                disc_clocks[sym] = i
                disc_clocks[ttsym] = i
                disc_clocks[rsym] = i
                disc_clocks[rttsym] = i

                insert_by_type!(disc_buffers[i], sym)
            end
            t = get_iv(sys)
            for eq in observed(disc_sys)
                # TODO: Is this a valid check
                # FIXME: This shouldn't be necessary
                eq.rhs === -0.0 && continue
                sym = eq.lhs
                ttsym = default_toterm(sym)
                rsym = renamespace(sys, sym)
                rttsym = renamespace(sys, ttsym)
                if iscall(sym) && operation(sym) == Shift(t, 1)
                    sym = only(arguments(sym))
                end
                disc_clocks[sym] = i
                disc_clocks[ttsym] = i
                disc_clocks[rsym] = i
                disc_clocks[rttsym] = i
            end
        end

        for par in inputs[continuous_id]
            is_parameter(sys, par) || error("Discrete subsystem input is not a parameter")
            par = unwrap(par)
            ttpar = default_toterm(par)
            rpar = renamespace(sys, par)
            rttpar = renamespace(sys, ttpar)
            iscall(par) && operation(par) isa Hold ||
                error("Continuous subsystem input is not a Hold")
            if haskey(disc_clocks, par)
                sym = par
            else
                sym = first(arguments(par))
            end
            haskey(disc_clocks, sym) ||
                error("Variable $par not part of a discrete subsystem")
            disc_clocks[par] = disc_clocks[sym]
            disc_clocks[ttpar] = disc_clocks[sym]
            disc_clocks[rpar] = disc_clocks[sym]
            disc_clocks[rttpar] = disc_clocks[sym]
            insert_by_type!(disc_buffers[disc_clocks[sym]], par)
        end
    end

    affs = vcat(affects(continuous_events(sys)), affects(discrete_events(sys)))
    user_affect_clock = maximum(values(disc_clocks); init = 0) + 1
    for affect in affs
        if affect isa Equation
            is_parameter(sys, affect.lhs) || continue
            sym = affect.lhs
            ttsym = default_toterm(sym)
            rsym = renamespace(sys, sym)
            rttsym = renamespace(sys, ttsym)

            disc_clocks[sym] = user_affect_clock
            disc_clocks[ttsym] = user_affect_clock
            disc_clocks[rsym] = user_affect_clock
            disc_clocks[rttsym] = user_affect_clock

            buffer = get!(disc_buffers, user_affect_clock, Dict{Any, Set{BasicSymbolic}}())
            insert_by_type!(buffer, affect.lhs)
        else
            discs = discretes(affect)
            for disc in discs
                is_parameter(sys, disc) ||
                    error("Expected discrete variable $disc in callback to be a parameter")
                disc = unwrap(disc)
                ttdisc = default_toterm(disc)
                rdisc = renamespace(sys, disc)
                rttdisc = renamespace(sys, ttdisc)
                disc_clocks[disc] = user_affect_clock
                disc_clocks[ttdisc] = user_affect_clock
                disc_clocks[rdisc] = user_affect_clock
                disc_clocks[rttdisc] = user_affect_clock

                buffer = get!(
                    disc_buffers, user_affect_clock, Dict{Any, Set{BasicSymbolic}}())
                insert_by_type!(buffer, disc)
            end
        end
    end

    if has_parameter_dependencies(sys)
        pdeps = parameter_dependencies(sys)
        for (sym, value) in pdeps
            sym = unwrap(sym)
            insert_by_type!(dependent_buffers, sym)
        end
    end

    for p in parameters(sys)
        p = unwrap(p)
        ctype = symtype(p)
        haskey(disc_clocks, p) && continue
        haskey(dependent_buffers, ctype) && p in dependent_buffers[ctype] && continue
        insert_by_type!(
            if ctype <: Real || ctype <: AbstractArray{<:Real}
                if istunable(p, true) && Symbolics.shape(p) !== Symbolics.Unknown()
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

    disc_idxs = Dict{Union{Symbol, BasicSymbolic}, Tuple{Int, Int, Int}}()
    disc_buffer_sizes = [BufferTemplate[] for _ in 1:length(disc_buffers)]
    disc_buffer_types = Set()
    for buffer in values(disc_buffers)
        union!(disc_buffer_types, keys(buffer))
    end

    for (clockidx, buffer) in disc_buffers
        for (i, btype) in enumerate(disc_buffer_types)
            if !haskey(buffer, btype)
                push!(disc_buffer_sizes[clockidx], BufferTemplate(btype, 0))
                continue
            end
            push!(disc_buffer_sizes[clockidx], BufferTemplate(btype, length(buffer[btype])))
            for (j, sym) in enumerate(buffer[btype])
                disc_idxs[sym] = (clockidx, i, j)
                disc_idxs[default_toterm(sym)] = (clockidx, i, j)
            end
        end
    end
    for (sym, clockid) in disc_clocks
        haskey(disc_idxs, sym) && continue
        disc_idxs[sym] = (clockid, 0, 0)
        disc_idxs[default_toterm(sym)] = (clockid, 0, 0)
    end

    function get_buffer_sizes_and_idxs(buffers::Dict{Any, Set{BasicSymbolic}})
        idxs = ParamIndexMap()
        buffer_sizes = BufferTemplate[]
        for (i, (T, buf)) in enumerate(buffers)
            for (j, p) in enumerate(buf)
                ttp = default_toterm(p)
                rp = renamespace(sys, p)
                rttp = renamespace(sys, ttp)
                idxs[p] = (i, j)
                idxs[ttp] = (i, j)
                idxs[rp] = (i, j)
                idxs[rttp] = (i, j)
            end
            push!(buffer_sizes, BufferTemplate(T, length(buf)))
        end
        return idxs, buffer_sizes
    end

    tunable_idxs, tunable_buffer_sizes = get_buffer_sizes_and_idxs(tunable_buffers)
    const_idxs, const_buffer_sizes = get_buffer_sizes_and_idxs(constant_buffers)
    dependent_idxs, dependent_buffer_sizes = get_buffer_sizes_and_idxs(dependent_buffers)
    nonnumeric_idxs, nonnumeric_buffer_sizes = get_buffer_sizes_and_idxs(nonnumeric_buffers)

    for sym in Iterators.flatten((keys(unk_idxs), keys(disc_idxs), keys(tunable_idxs),
        keys(const_idxs), keys(dependent_idxs), keys(nonnumeric_idxs),
        observed_syms, independent_variable_symbols(sys)))
        if hasname(sym) && (!iscall(sym) || operation(sym) !== getindex)
            symbol_to_variable[getname(sym)] = sym
        end
    end

    return IndexCache(
        unk_idxs,
        disc_idxs,
        tunable_idxs,
        const_idxs,
        dependent_idxs,
        nonnumeric_idxs,
        observed_syms,
        disc_buffer_sizes,
        tunable_buffer_sizes,
        const_buffer_sizes,
        dependent_buffer_sizes,
        nonnumeric_buffer_sizes,
        symbol_to_variable
    )
end

function SymbolicIndexingInterface.is_variable(ic::IndexCache, sym)
    if sym isa Symbol
        sym = get(ic.symbol_to_variable, sym, nothing)
        sym === nothing && return false
    end
    return check_index_map(ic.unknown_idx, sym) !== nothing
end

function SymbolicIndexingInterface.variable_index(ic::IndexCache, sym)
    if sym isa Symbol
        sym = get(ic.symbol_to_variable, sym, nothing)
        sym === nothing && return nothing
    end
    return check_index_map(ic.unknown_idx, sym)
end

function SymbolicIndexingInterface.is_parameter(ic::IndexCache, sym)
    if sym isa Symbol
        sym = get(ic.symbol_to_variable, sym, nothing)
        sym === nothing && return false
    end
    return check_index_map(ic.tunable_idx, sym) !== nothing ||
           check_index_map(ic.discrete_idx, sym) !== nothing ||
           check_index_map(ic.constant_idx, sym) !== nothing ||
           check_index_map(ic.nonnumeric_idx, sym) !== nothing ||
           check_index_map(ic.dependent_idx, sym) !== nothing
end

function SymbolicIndexingInterface.parameter_index(ic::IndexCache, sym)
    if sym isa Symbol
        sym = get(ic.symbol_to_variable, sym, nothing)
        sym === nothing && return nothing
    end
    validate_size = Symbolics.isarraysymbolic(sym) &&
                    Symbolics.shape(sym) !== Symbolics.Unknown()
    return if (idx = check_index_map(ic.tunable_idx, sym)) !== nothing
        ParameterIndex(SciMLStructures.Tunable(), idx, validate_size)
    elseif (idx = check_index_map(ic.discrete_idx, sym)) !== nothing
        ParameterIndex(SciMLStructures.Discrete(), idx, validate_size)
    elseif (idx = check_index_map(ic.constant_idx, sym)) !== nothing
        ParameterIndex(SciMLStructures.Constants(), idx, validate_size)
    elseif (idx = check_index_map(ic.nonnumeric_idx, sym)) !== nothing
        ParameterIndex(NONNUMERIC_PORTION, idx, validate_size)
    elseif (idx = check_index_map(ic.dependent_idx, sym)) !== nothing
        ParameterIndex(DEPENDENT_PORTION, idx, validate_size)
    else
        nothing
    end
end

function SymbolicIndexingInterface.is_timeseries_parameter(ic::IndexCache, sym)
    if sym isa Symbol
        sym = get(ic.symbol_to_variable, sym, nothing)
        sym === nothing && return false
    end
    return check_index_map(ic.discrete_idx, sym) !== nothing
end

function SymbolicIndexingInterface.timeseries_parameter_index(ic::IndexCache, sym)
    if sym isa Symbol
        sym = get(ic.symbol_to_variable, sym, nothing)
        sym === nothing && return nothing
    end
    idx = check_index_map(ic.discrete_idx, sym)
    idx === nothing && return nothing
    clockid, partitionid... = idx
    return ParameterTimeseriesIndex(clockid, partitionid)
end

function check_index_map(idxmap, sym)
    if (idx = get(idxmap, sym, nothing)) !== nothing
        return idx
    elseif !isa(sym, Symbol) && (!iscall(sym) || operation(sym) !== getindex) &&
           hasname(sym) && (idx = get(idxmap, getname(sym), nothing)) !== nothing
        return idx
    end
    dsym = default_toterm(sym)
    isequal(sym, dsym) && return nothing
    if (idx = get(idxmap, dsym, nothing)) !== nothing
        idx
    elseif !isa(dsym, Symbol) && (!iscall(dsym) || operation(dsym) !== getindex) &&
           hasname(dsym) && (idx = get(idxmap, getname(dsym), nothing)) !== nothing
        idx
    else
        nothing
    end
end

function discrete_linear_index(ic::IndexCache, idx::ParameterIndex)
    idx.portion isa SciMLStructures.Discrete || error("Discrete variable index expected")
    ind = sum(temp.length for temp in ic.tunable_buffer_sizes; init = 0)
    for clockbuftemps in Iterators.take(ic.discrete_buffer_sizes, idx.idx[1] - 1)
        ind += sum(temp.length for temp in clockbuftemps; init = 0)
    end
    ind += sum(
        temp.length
        for temp in Iterators.take(ic.discrete_buffer_sizes[idx.idx[1]], idx.idx[2] - 1);
        init = 0)
    ind += idx.idx[3]
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
    for temp in Iterators.flatten(ic.discrete_buffer_sizes))
    const_buf = Tuple(BasicSymbolic[unwrap(variable(:DEF)) for _ in 1:(temp.length)]
    for temp in ic.constant_buffer_sizes)
    dep_buf = Tuple(BasicSymbolic[unwrap(variable(:DEF)) for _ in 1:(temp.length)]
    for temp in ic.dependent_buffer_sizes)
    nonnumeric_buf = Tuple(BasicSymbolic[unwrap(variable(:DEF)) for _ in 1:(temp.length)]
    for temp in ic.nonnumeric_buffer_sizes)
    for p in ps
        p = unwrap(p)
        if haskey(ic.discrete_idx, p)
            disc_offset = length(first(ic.discrete_buffer_sizes))
            i, j, k = ic.discrete_idx[p]
            disc_buf[(i - 1) * disc_offset + j][k] = p
        elseif haskey(ic.tunable_idx, p)
            i, j = ic.tunable_idx[p]
            param_buf[i][j] = p
        elseif haskey(ic.constant_idx, p)
            i, j = ic.constant_idx[p]
            const_buf[i][j] = p
        elseif haskey(ic.dependent_idx, p)
            i, j = ic.dependent_idx[p]
            dep_buf[i][j] = p
        elseif haskey(ic.nonnumeric_idx, p)
            i, j = ic.nonnumeric_idx[p]
            nonnumeric_buf[i][j] = p
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
