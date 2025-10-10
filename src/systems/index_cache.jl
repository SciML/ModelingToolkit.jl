const TypeT = Union{DataType, UnionAll, Union}

struct BufferTemplate
    type::TypeT
    length::Int
end

struct Nonnumeric <: SciMLStructures.AbstractPortion end
const NONNUMERIC_PORTION = Nonnumeric()

struct ParameterIndex{P, I}
    portion::P
    idx::I
    validate_size::Bool
end

ParameterIndex(portion, idx) = ParameterIndex(portion, idx, false)
ParameterIndex(p::ParameterIndex) = ParameterIndex(p.portion, p.idx, false)

struct DiscreteIndex
    # of all buffers corresponding to types, which one
    buffer_idx::Int
    # Index in the above buffer
    idx_in_buffer::Int
    # Which clock (corresponds to Block of BlockedArray)
    clock_idx::Int
    # Which index in `buffer[Block(clockidx)]`
    idx_in_clock::Int
end

const MaybeUnknownArrayIndexT = Union{Int, UnitRange{Int}, AbstractArray{Int}}
const MaybeArrayIndexT = Union{Int, UnitRange{Int}, Base.ReshapedArray{Int, N, UnitRange{Int}} where {N}}
const ParamIndexMap = Dict{SymbolicT, Tuple{Int, Int}}
const NonnumericMap = Dict{SymbolicT, Tuple{Int, Int}}
const UnknownIndexMap = Dict{SymbolicT, MaybeUnknownArrayIndexT}
const TunableIndexMap = Dict{SymbolicT, MaybeArrayIndexT}
const TimeseriesSetType = Set{Union{ContinuousTimeseries, Int}}

const SymbolicParam = SymbolicT

struct IndexCache
    unknown_idx::UnknownIndexMap
    # sym => (bufferidx, idx_in_buffer)
    discrete_idx::Dict{SymbolicParam, DiscreteIndex}
    # sym => (clockidx, idx_in_clockbuffer)
    callback_to_clocks::Dict{Any, Vector{Int}}
    tunable_idx::TunableIndexMap
    initials_idx::TunableIndexMap
    constant_idx::ParamIndexMap
    nonnumeric_idx::NonnumericMap
    observed_syms_to_timeseries::Dict{SymbolicT, TimeseriesSetType}
    dependent_pars_to_timeseries::Dict{SymbolicT, TimeseriesSetType}
    discrete_buffer_sizes::Vector{Vector{BufferTemplate}}
    tunable_buffer_size::BufferTemplate
    initials_buffer_size::BufferTemplate
    constant_buffer_sizes::Vector{BufferTemplate}
    nonnumeric_buffer_sizes::Vector{BufferTemplate}
    symbol_to_variable::Dict{Symbol, SymbolicParam}
end

function IndexCache(sys::AbstractSystem)
    unks = unknowns(sys)
    unk_idxs = UnknownIndexMap()
    symbol_to_variable = Dict{Symbol, SymbolicParam}()

    let idx = 1
        for sym in unks
            rsym = renamespace(sys, sym)
            sym_idx::MaybeUnknownArrayIndexT = if Symbolics.isarraysymbolic(sym)
                reshape(idx:(idx + length(sym) - 1), size(sym))
            else
                idx
            end
            unk_idxs[sym] = sym_idx
            unk_idxs[rsym] = sym_idx
            idx += length(sym)
        end
        found_array_syms = Set{SymbolicT}()
        for sym in unks
            iscall(sym) && operation(sym) === getindex || continue
            arrsym = arguments(sym)[1]
            arrsym in found_array_syms && continue
            idxs = Int[]
            valid_arrsym = true
            for i in eachindex(arrsym)
                idxsym = arrsym[i]
                idx = get(unk_idxs, idxsym, nothing)::Union{Int, Nothing}
                valid_arrsym = idx !== nothing
                valid_arrsym || break
                push!(idxs, idx::Int)
            end
            push!(found_array_syms, arrsym)
            valid_arrsym || break
            if idxs == idxs[begin]:idxs[end]
                idxs = reshape(idxs[begin]:idxs[end], size(idxs))::AbstractArray{Int}
            else
                idxs = reshape(idxs, size(arrsym))::AbstractArray{Int}
            end
            rsym = renamespace(sys, arrsym)
            unk_idxs[arrsym] = idxs
            unk_idxs[rsym] = idxs
        end
    end

    tunable_pars = SymbolicT[]
    initial_pars = SymbolicT[]
    constant_buffers = Dict{TypeT, Set{SymbolicT}}()
    nonnumeric_buffers = Dict{TypeT, Set{SymbolicT}}()

    disc_param_callbacks = Dict{SymbolicParam, BitSet}()
    cevs = continuous_events(sys)
    devs = discrete_events(sys)
    events = Union{SymbolicContinuousCallback, SymbolicDiscreteCallback}[cevs; devs]
    parse_callbacks_for_discretes!(cevs, disc_param_callbacks, constant_buffers, nonnumeric_buffers, 0)
    parse_callbacks_for_discretes!(devs, disc_param_callbacks, constant_buffers, nonnumeric_buffers, length(cevs))
    clock_partitions = unique(collect(values(disc_param_callbacks)))::Vector{BitSet}
    disc_symtypes = Set{TypeT}()
    for x in keys(disc_param_callbacks)
        push!(disc_symtypes, symtype(x))
    end
    disc_symtypes = collect(disc_symtypes)::Vector{TypeT}
    disc_symtype_idx = Dict{TypeT, Int}(zip(disc_symtypes, eachindex(disc_symtypes)))
    disc_syms_by_symtype = [SymbolicParam[] for _ in disc_symtypes]
    for sym in keys(disc_param_callbacks)
        push!(disc_syms_by_symtype[disc_symtype_idx[symtype(sym)]], sym)
    end
    disc_syms_by_symtype_by_partition = [Vector{SymbolicParam}[] for _ in disc_symtypes]
    for (i, buffer) in enumerate(disc_syms_by_symtype)
        for partition in clock_partitions
            push!(disc_syms_by_symtype_by_partition[i], filter(==(partition) ∘ Base.Fix1(getindex, disc_param_callbacks), buffer))
        end
    end
    disc_idxs = Dict{SymbolicParam, DiscreteIndex}()
    callback_to_clocks = Dict{
        Union{SymbolicContinuousCallback, SymbolicDiscreteCallback}, BitSet}()
    for (typei, disc_syms_by_partition) in enumerate(disc_syms_by_symtype_by_partition)
        symi = 0
        for (parti, disc_syms) in enumerate(disc_syms_by_partition)
            for clockidx in clock_partitions[parti]
                buffer = get!(() -> Set{Int}(), callback_to_clocks, events[clockidx])
                push!(buffer, parti)
            end
            clocki = 0
            for sym in disc_syms
                symi += 1
                clocki += 1
                ttsym = default_toterm(sym)
                rsym = renamespace(sys, sym)
                rttsym = renamespace(sys, ttsym)
                for cursym in (sym, ttsym, rsym, rttsym)
                    disc_idxs[cursym] = DiscreteIndex(typei, symi, parti, clocki)
                end
            end
        end
    end
    callback_to_clocks = Dict{
        Union{SymbolicContinuousCallback, SymbolicDiscreteCallback}, Vector{Int}}(k => collect(v)
    for (k, v) in callback_to_clocks)

    disc_buffer_templates = Vector{BufferTemplate}[]
    for (symtype, disc_syms_by_partition) in zip(
        disc_symtypes, disc_syms_by_symtype_by_partition)
        push!(disc_buffer_templates, map(Base.Fix1(BufferTemplate, symtype) ∘ length, disc_syms_by_partition))
    end

    for p in parameters(sys; initial_parameters = true)
        ctype = symtype(p)
        if ctype <: FnType
            ctype = fntype_to_function_type(ctype)::TypeT
        end
        haskey(disc_idxs, p) && continue
        haskey(constant_buffers, ctype) && p in constant_buffers[ctype] && continue
        haskey(nonnumeric_buffers, ctype) && p in nonnumeric_buffers[ctype] && continue
        insert_by_type!(
            if ctype <: Real || ctype <: AbstractArray{<:Real}
                if istunable(p, true) && symbolic_has_known_size(p) &&
                   (ctype == Real || ctype <: AbstractFloat ||
                    ctype <: AbstractArray{Real} ||
                    ctype <: AbstractArray{<:AbstractFloat})
                    if iscall(p) && operation(p) === Initial()
                        initial_pars
                    else
                        tunable_pars
                    end
                else
                    constant_buffers
                end
            else
                nonnumeric_buffers
            end,
            p,
            ctype
        )
    end

    const_idxs,
    const_buffer_sizes = get_buffer_sizes_and_idxs(ParamIndexMap, sys, constant_buffers)
    nonnumeric_idxs,
    nonnumeric_buffer_sizes = get_buffer_sizes_and_idxs(NonnumericMap, sys, nonnumeric_buffers)

    tunable_idxs = TunableIndexMap()
    tunable_buffer_size = 0
    if is_initializesystem(sys)
        append!(tunable_pars, initial_pars)
        empty!(initial_pars)
    end
    for p in tunable_pars
        sh = SU.shape(p)
        idx = if !SU.is_array_shape(sh)
            tunable_buffer_size + 1
        else
            reshape(
                (tunable_buffer_size + 1):(tunable_buffer_size + length(p)), size(p))
        end
        tunable_buffer_size += length(p)
        tunable_idxs[p] = idx
        tunable_idxs[default_toterm(p)] = idx
        if hasname(p) && (!iscall(p) || operation(p) !== getindex)
            symbol_to_variable[getname(p)] = p
            symbol_to_variable[getname(default_toterm(p))] = p
        end
    end

    initials_idxs = TunableIndexMap()
    initials_buffer_size = 0
    for p in initial_pars
        sh = SU.shape(p)
        idx = if !SU.is_array_shape(sh)
            initials_buffer_size + 1
        else
            reshape(
                (initials_buffer_size + 1):(initials_buffer_size + length(p)), size(p))
        end
        initials_buffer_size += length(p)
        initials_idxs[p] = idx
        initials_idxs[default_toterm(p)] = idx
        if hasname(p) && (!iscall(p) || operation(p) !== getindex)
            symbol_to_variable[getname(p)] = p
            symbol_to_variable[getname(default_toterm(p))] = p
        end
    end

    for k in collect(keys(tunable_idxs))
        v = tunable_idxs[k]
        v isa AbstractArray || continue
        v = v::Union{UnitRange{Int}, Base.ReshapedArray{Int, N, UnitRange{Int}} where {N}}
        iter = vec(collect(k)::Array{SymbolicT})::Vector{SymbolicT}
        for (kk::SymbolicT, vv) in zip(iter, v)
            tunable_idxs[kk] = vv
        end
    end
    for k in collect(keys(initials_idxs))
        v = initials_idxs[k]
        v isa AbstractArray || continue
        v = v::Union{UnitRange{Int}, Base.ReshapedArray{Int, N, UnitRange{Int}} where {N}}
        iter = vec(collect(k)::Array{SymbolicT})::Vector{SymbolicT}
        for (kk, vv) in zip(iter, v)
            initials_idxs[kk] = vv
        end
    end

    dependent_pars_to_timeseries = Dict{SymbolicT, TimeseriesSetType}()
    vs = Set{SymbolicT}()
    for eq in get_parameter_dependencies(sys)
        sym = eq.lhs
        SU.search_variables!(vs, eq.rhs)
        timeseries = TimeseriesSetType()
        if is_time_dependent(sys)
            for v in vs
                if (idx = get(disc_idxs, v, nothing)) !== nothing
                    push!(timeseries, idx.clock_idx)
                end
            end
        end
        ttsym = default_toterm(sym)
        rsym = renamespace(sys, sym)
        rttsym = renamespace(sys, ttsym)
        for s in (sym, ttsym, rsym, rttsym)
            dependent_pars_to_timeseries[s] = timeseries
            if hasname(s) && (!iscall(s) || operation(s) !== getindex)
                symbol_to_variable[getname(s)] = sym
            end
        end
    end

    observed_syms_to_timeseries = Dict{SymbolicT, TimeseriesSetType}()
    for eq in observed(sys)
        if symbolic_type(eq.lhs) != NotSymbolic()
            sym = eq.lhs
            empty!(vs)
            SU.search_variables!(vs, eq.rhs)
            timeseries = TimeseriesSetType()
            if is_time_dependent(sys)
                for v in vs
                    if (idx = get(disc_idxs, v, nothing)) !== nothing
                        push!(timeseries, idx.clock_idx)
                    elseif Moshi.Match.@match v begin
                            BSImpl.Term(; f, args) => begin
                                f === getindex && (idx = get(disc_idxs, args[1], nothing)) !== nothing
                            end
                            _ => false
                        end
                        push!(timeseries, idx.clock_idx)
                    elseif haskey(observed_syms_to_timeseries, v)
                        union!(timeseries, observed_syms_to_timeseries[v])
                    elseif haskey(dependent_pars_to_timeseries, v)
                        union!(timeseries, dependent_pars_to_timeseries[v])
                    end
                end
                if isempty(timeseries)
                    push!(timeseries, ContinuousTimeseries())
                end
            end
            ttsym = default_toterm(sym)
            rsym = renamespace(sys, sym)
            rttsym = renamespace(sys, ttsym)
            for s in (sym, ttsym, rsym, rttsym)
                observed_syms_to_timeseries[s] = timeseries
            end
        end
    end

    populate_symbol_to_var!(symbol_to_variable, keys(unk_idxs))
    populate_symbol_to_var!(symbol_to_variable, keys(disc_idxs))
    populate_symbol_to_var!(symbol_to_variable, keys(tunable_idxs))
    populate_symbol_to_var!(symbol_to_variable, keys(const_idxs))
    populate_symbol_to_var!(symbol_to_variable, keys(nonnumeric_idxs))
    populate_symbol_to_var!(symbol_to_variable, independent_variable_symbols(sys))

    return IndexCache(
        unk_idxs,
        disc_idxs,
        callback_to_clocks,
        tunable_idxs,
        initials_idxs,
        const_idxs,
        nonnumeric_idxs,
        observed_syms_to_timeseries,
        dependent_pars_to_timeseries,
        disc_buffer_templates,
        BufferTemplate(Number, tunable_buffer_size),
        BufferTemplate(Number, initials_buffer_size),
        const_buffer_sizes,
        nonnumeric_buffer_sizes,
        symbol_to_variable
    )
end

function populate_symbol_to_var!(symbol_to_variable::Dict{Symbol, SymbolicT}, vars)
    for sym::SymbolicT in vars
        if hasname(sym) && (!iscall(sym) || operation(sym) !== getindex)
            symbol_to_variable[getname(sym)] = sym
        end
    end
end

"""
    $TYPEDSIGNATURES

Utility function for the `IndexCache` constructor.
"""
function insert_by_type!(buffers::Dict{TypeT, Set{SymbolicT}}, sym::SymbolicT, ctype::TypeT)
    buf = get!(Set{SymbolicT}, buffers, ctype)
    push!(buf, sym)
end
function insert_by_type!(buffers::Vector{SymbolicT}, sym::SymbolicT, ::TypeT)
    push!(buffers, sym)
end

function parse_callbacks_for_discretes!(events::Vector, disc_param_callbacks::Dict{SymbolicT, BitSet}, constant_buffers::Dict{TypeT, Set{SymbolicT}}, nonnumeric_buffers::Dict{TypeT, Set{SymbolicT}}, offset::Int)
    for (i, event) in enumerate(events)
        discs = Set{SymbolicParam}()
        affect = event.affect::Union{AffectSystem, ImperativeAffect, Nothing}
        if affect isa AffectSystem || affect isa ImperativeAffect
            union!(discs, discretes(affect))
        elseif affect === nothing
            continue
        end

        for sym in discs
            is_parameter(sys, sym) ||
                error("Expected discrete variable $sym in callback to be a parameter")

            # Only `foo(t)`-esque parameters can be saved
            if iscall(sym) && length(arguments(sym)) == 1 &&
               isequal(only(arguments(sym)), get_iv(sys))
                clocks = get!(BitSet, disc_param_callbacks, sym)
                push!(clocks, i + offset)
            elseif is_variable_floatingpoint(sym)
                insert_by_type!(constant_buffers, sym, symtype(sym))
            else
                stype = symtype(sym)
                if stype <: FnType
                    stype = fntype_to_function_type(stype)::TypeT
                end
                insert_by_type!(nonnumeric_buffers, sym, stype)
            end
        end
    end
end

function get_buffer_sizes_and_idxs(::Type{BufT}, sys::AbstractSystem, buffers::Dict) where {BufT}
    idxs = BufT()
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
        if T <: Symbolics.FnType
            T = Any
        end
        push!(buffer_sizes, BufferTemplate(T, length(buf)))
    end
    return idxs, buffer_sizes
end

function SymbolicIndexingInterface.is_variable(ic::IndexCache, sym)
    variable_index(ic, sym) !== nothing
end

function SymbolicIndexingInterface.variable_index(ic::IndexCache, sym::Union{Num, Symbolics.Arr, Symbolics.CallAndWrap})
    variable_index(ic, unwrap(sym))
end
function SymbolicIndexingInterface.variable_index(ic::IndexCache, sym::Symbol)
    sym = get(ic.symbol_to_variable, sym, nothing)
    sym === nothing && return nothing
    variable_index(ic, sym)
end
function SymbolicIndexingInterface.variable_index(ic::IndexCache, sym::SymbolicT)
    idx = check_index_map(ic.unknown_idx, sym)
    idx === nothing || return idx
    iscall(sym) && operation(sym) == getindex || return nothing
    args = arguments(sym)
    idx = variable_index(ic, args[1])
    idx === nothing && return nothing
    return idx[args[2:end]...]
end
SymbolicIndexingInterface.variable_index(ic::IndexCache, sym) = false

function SymbolicIndexingInterface.is_parameter(ic::IndexCache, sym)
    parameter_index(ic, sym) !== nothing
end

function SymbolicIndexingInterface.parameter_index(ic::IndexCache, sym::Union{Num, Symbolics.Arr, Symbolics.CallAndWrap})
    parameter_index(ic, unwrap(sym))
end
function SymbolicIndexingInterface.parameter_index(ic::IndexCache, sym::Symbol)
    sym = get(ic.symbol_to_variable, sym, nothing)
    sym === nothing && return nothing
    parameter_index(ic, sym)
end
function SymbolicIndexingInterface.parameter_index(ic::IndexCache, sym::SymbolicT)
    validate_size = Symbolics.isarraysymbolic(sym) && symbolic_has_known_size(sym)
    return if (idx = check_index_map(ic.tunable_idx, sym)) !== nothing
        ParameterIndex(SciMLStructures.Tunable(), idx, validate_size)
    elseif (idx = check_index_map(ic.initials_idx, sym)) !== nothing
        ParameterIndex(SciMLStructures.Initials(), idx, validate_size)
    elseif (idx = check_index_map(ic.discrete_idx, sym)) !== nothing
        ParameterIndex(
            SciMLStructures.Discrete(), (idx.buffer_idx, idx.idx_in_buffer), validate_size)
    elseif (idx = check_index_map(ic.constant_idx, sym)) !== nothing
        ParameterIndex(SciMLStructures.Constants(), idx, validate_size)
    elseif (idx = check_index_map(ic.nonnumeric_idx, sym)) !== nothing
        ParameterIndex(NONNUMERIC_PORTION, idx, validate_size)
    elseif iscall(sym) && operation(sym) == getindex
        args = arguments(sym)
        pidx = parameter_index(ic, args[1])
        pidx === nothing && return nothing
        if pidx.portion == SciMLStructures.Tunable()
            ParameterIndex(pidx.portion,
                Origin(first.(axes((args[1]))))(reshape(pidx.idx, size(args[1])))[args[2:end]...],
                pidx.validate_size)
        else
            ParameterIndex(pidx.portion, (pidx.idx..., args[2:end]...), pidx.validate_size)
        end
    end
end

function SymbolicIndexingInterface.is_timeseries_parameter(ic::IndexCache, sym)
    timeseries_parameter_index(ic, sym) !== nothing
end

function SymbolicIndexingInterface.timeseries_parameter_index(ic::IndexCache, sym)
    if sym isa Symbol
        sym = get(ic.symbol_to_variable, sym, nothing)
        sym === nothing && return nothing
    end
    sym = unwrap(sym)
    idx = check_index_map(ic.discrete_idx, sym)
    idx === nothing ||
        return ParameterTimeseriesIndex(idx.clock_idx, (idx.buffer_idx, idx.idx_in_clock))
    iscall(sym) && operation(sym) == getindex || return nothing
    args = arguments(sym)
    idx = timeseries_parameter_index(ic, args[1])
    idx === nothing && return nothing
    return ParameterTimeseriesIndex(
        idx.timeseries_idx, (idx.parameter_idx..., args[2:end]...))
end

function check_index_map(idxmap::Dict{SymbolicT, V}, sym::SymbolicT)::Union{V, Nothing} where {V}
    idx = get(idxmap, sym, nothing)
    idx === nothing || return idx
    dsym = default_toterm(sym)
    isequal(sym, dsym) && return nothing
    idx = get(idxmap, dsym, nothing)
    idx === nothing || return idx
    return nothing
end

const ReorderedParametersT = Vector{Union{Vector{SymbolicT}, Vector{Vector{SymbolicT}}}}

function reorder_parameters(
        sys::AbstractSystem, ps = parameters(sys; initial_parameters = true); kwargs...)
    if has_index_cache(sys) && get_index_cache(sys) !== nothing
        return reorder_parameters(get_index_cache(sys)::IndexCache, ps; kwargs...)
    elseif ps isa Tuple
        return ReorderedParametersT(collect(ps))
    else
        return eltype(ReorderedParametersT)[ps]
    end
end

const COMMON_DEFAULT_VAR = unwrap(only(@variables __DEF__))

function reorder_parameters(ic::IndexCache, ps::Vector{SymbolicT}; drop_missing = false, flatten = true)
    result = ReorderedParametersT()
    isempty(ps) && return result
    param_buf = fill(COMMON_DEFAULT_VAR, ic.tunable_buffer_size.length)
    if !isempty(param_buf) || !flatten
        push!(result, param_buf)
    end
    initials_buf = fill(COMMON_DEFAULT_VAR, ic.initials_buffer_size.length)
    if !isempty(initials_buf) || !flatten
        push!(result, initials_buf)
    end

    disc_buf = Vector{SymbolicT}[]
    for bufszs in ic.discrete_buffer_sizes
        push!(disc_buf, fill(COMMON_DEFAULT_VAR, sum(x -> x.length, bufszs)))
    end
    const_buf = Vector{SymbolicT}[]
    for bufsz in ic.constant_buffer_sizes
        push!(const_buf, fill(COMMON_DEFAULT_VAR, bufsz.length))
    end
    nonnumeric_buf = Vector{SymbolicT}[]
    for bufsz in ic.nonnumeric_buffer_sizes
        push!(nonnumeric_buf, fill(COMMON_DEFAULT_VAR, bufsz.length))
    end
    if flatten
        append!(result, disc_buf)
        append!(result, const_buf)
        append!(result, nonnumeric_buf)
    else
        push!(result, disc_buf)
        push!(result, const_buf)
        push!(result, nonnumeric_buf)
    end
    for p in ps
        if haskey(ic.discrete_idx, p)
            idx = ic.discrete_idx[p]
            disc_buf[idx.buffer_idx][idx.idx_in_buffer] = p
        elseif haskey(ic.tunable_idx, p)
            i = ic.tunable_idx[p]
            if i isa Int
                param_buf[i] = p
            else
                param_buf[i] = collect(p)
            end
        elseif haskey(ic.initials_idx, p)
            i = ic.initials_idx[p]
            if i isa Int
                initials_buf[i] = p
            else
                initials_buf[i] = collect(p)
            end
        elseif haskey(ic.constant_idx, p)
            i, j = ic.constant_idx[p]
            const_buf[i][j] = p
        elseif haskey(ic.nonnumeric_idx, p)
            i, j = ic.nonnumeric_idx[p]
            nonnumeric_buf[i][j] = p
        else
            error("Invalid parameter $p")
        end
    end

    if drop_missing
        filterer = !isequal(COMMON_DEFAULT_VAR)
        for inner in result
            if inner isa Vector{SymbolicT}
                filter!(filterer, inner)
            elseif inner isa Vector{Vector{SymbolicT}}
                for buf in inner
                    filter!(filterer, buf)
                end
            end
        end
    end

    return result
end

# Given a parameter index, find the index of the buffer it is in when
# `MTKParameters` is iterated
function iterated_buffer_index(ic::IndexCache, ind::ParameterIndex)
    idx = 0
    if ind.portion isa SciMLStructures.Tunable
        return idx + 1
    elseif ic.tunable_buffer_size.length > 0
        idx += 1
    end
    if ind.portion isa SciMLStructures.Initials
        return idx + 1
    elseif ic.initials_buffer_size.length > 0
        idx += 1
    end
    if ind.portion isa SciMLStructures.Discrete
        return idx + ind.idx[1]
    elseif !isempty(ic.discrete_buffer_sizes)
        idx += length(ic.discrete_buffer_sizes)
    end
    if ind.portion isa SciMLStructures.Constants
        return idx + ind.idx[1]
    elseif !isempty(ic.constant_buffer_sizes)
        idx += length(ic.constant_buffer_sizes)
    end
    if ind.portion == NONNUMERIC_PORTION
        return idx + ind.idx[1]
    end
    error("Unhandled portion $(ind.portion)")
end

function get_buffer_template(ic::IndexCache, pidx::ParameterIndex)
    (; portion, idx) = pidx

    if portion isa SciMLStructures.Tunable
        return ic.tunable_buffer_size
    elseif portion isa SciMLStructures.Initials
        return ic.initials_buffer_size
    elseif portion isa SciMLStructures.Discrete
        return ic.discrete_buffer_sizes[idx[1]][1]
    elseif portion isa SciMLStructures.Constants
        return ic.constant_buffer_sizes[idx[1]]
    elseif portion isa Nonnumeric
        return ic.nonnumeric_buffer_sizes[idx[1]]
    else
        error("Unhandled portion $portion")
    end
end

fntype_to_function_type(::Type{FnType{A, R, T}}) where {A, R, T} = T
fntype_to_function_type(::Type{FnType{A, R, Nothing}}) where {A, R} = FunctionWrapper{R, A}
fntype_to_function_type(::Type{FnType{A, R}}) where {A, R} = FunctionWrapper{R, A}

"""
    reorder_dimension_by_tunables!(dest::AbstractArray, sys::AbstractSystem, arr::AbstractArray, syms; dim = 1)

Assuming the order of values in dimension `dim` of `arr` correspond to the order of tunable
parameters in the system, reorder them according to the order described in `syms`. `syms` must
be a permutation of `tunable_parameters(sys)`. The result is written to `dest`. The `size` of `dest` and
`arr` must be equal. Return `dest`.

See also: [`MTKParameters`](@ref), [`tunable_parameters`](@ref), [`reorder_dimension_by_tunables`](@ref).
"""
function reorder_dimension_by_tunables!(
        dest::AbstractArray, sys::AbstractSystem, arr::AbstractArray, syms; dim = 1)
    if !iscomplete(sys)
        throw(ArgumentError("A completed system is required. Call `complete` or `mtkcompile` on the system."))
    end
    if !has_index_cache(sys) || (ic = get_index_cache(sys)) === nothing
        throw(ArgumentError("The system does not have an index cache. Call `complete` or `mtkcompile` on the system with `split = true`."))
    end
    if size(dest) != size(arr)
        throw(ArgumentError("Source and destination arrays must have the same size. Got source array with size $(size(arr)) and destination with size $(size(dest))."))
    end

    dsti = 1
    for sym in syms
        idx = parameter_index(ic, sym)
        if !(idx.portion isa SciMLStructures.Tunable)
            throw(ArgumentError("`syms` must be a permutation of `tunable_parameters(sys)`. Found $sym which is not a tunable parameter."))
        end

        dstidx = ntuple(
            i -> i == dim ? (dsti:(dsti + length(sym) - 1)) : (:), Val(ndims(arr)))
        destv = @view dest[dstidx...]
        dsti += length(sym)
        arridx = ntuple(i -> i == dim ? (idx.idx) : (:), Val(ndims(arr)))
        srcv = @view arr[arridx...]
        copyto!(destv, srcv)
    end
    return dest
end

"""
    reorder_dimension_by_tunables(sys::AbstractSystem, arr::AbstractArray, syms; dim = 1)

Out-of-place version of [`reorder_dimension_by_tunables!`](@ref).
"""
function reorder_dimension_by_tunables(
        sys::AbstractSystem, arr::AbstractArray, syms; dim = 1)
    buffer = similar(arr)
    reorder_dimension_by_tunables!(buffer, sys, arr, syms; dim)
    return buffer
end

function subset_unknowns_observed(
        ic::IndexCache, sys::AbstractSystem, newunknowns, newobsvars)
    unknown_idx = copy(ic.unknown_idx)
    empty!(unknown_idx)
    for (i, sym) in enumerate(newunknowns)
        ttsym = default_toterm(sym)
        rsym = renamespace(sys, sym)
        rttsym = renamespace(sys, ttsym)
        unknown_idx[sym] = unknown_idx[ttsym] = unknown_idx[rsym] = unknown_idx[rttsym] = i
    end
    observed_syms_to_timeseries = copy(ic.observed_syms_to_timeseries)
    empty!(observed_syms_to_timeseries)
    for sym in newobsvars
        ttsym = default_toterm(sym)
        rsym = renamespace(sys, sym)
        rttsym = renamespace(sys, ttsym)
        for s in (sym, ttsym, rsym, rttsym)
            observed_syms_to_timeseries[s] = ic.observed_syms_to_timeseries[sym]
        end
    end
    ic = @set ic.unknown_idx = unknown_idx
    @set! ic.observed_syms_to_timeseries = observed_syms_to_timeseries
    return ic
end
