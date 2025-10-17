const TypeT = Union{DataType, UnionAll, Union}
struct BufferTemplate
end
struct ParameterIndex{P, I}
end
struct DiscreteIndex
end
struct IndexCache
end
function IndexCache(sys::AbstractSystem)
    let idx = 1
        for sym in unks
            sym_idx::MaybeUnknownArrayIndexT = if Symbolics.isarraysymbolic(sym)
            end
        end
        for sym in unks
        end
    end
    for eq in get_parameter_dependencies(sys)
        if is_time_dependent(sys)
            for v in vs
                if (idx = get(disc_idxs, v, nothing)) !== nothing
                end
            end
        end
        for s in (sym, ttsym, rsym, rttsym)
            if hasname(s) && (!iscall(s) || operation(s) !== getindex)
            end
        end
    end
    for eq in observed(sys)
        if symbolic_type(eq.lhs) != NotSymbolic()
            if is_time_dependent(sys)
                for v in vs
                    if (idx = get(disc_idxs, v, nothing)) !== nothing
                    elseif Moshi.Match.@match v begin
                            BSImpl.Term(; f, args) => begin
                            end
                        end
                    end
                end
            end
        end
    end
    return IndexCache(
    )
    for sym::SymbolicT in vars
        if hasname(sym) && (!iscall(sym) || operation(sym) !== getindex)
        end
    end
end
function insert_by_type!(buffers::Vector{SymbolicT}, sym::SymbolicT, ::TypeT)
    for (i, event) in enumerate(events)
        if affect isa AffectSystem || affect isa ImperativeAffect
            if iscall(sym) && length(arguments(sym)) == 1 &&
                if stype <: FnType
                end
            end
        end
    end
end
function get_buffer_sizes_and_idxs(::Type{BufT}, sys::AbstractSystem, buffers::Dict) where {BufT}
    for (i, (T, buf)) in enumerate(buffers)
        for (j, p) in enumerate(buf)
        end
    end
    if !isempty(initials_buf) || !flatten
    end
    for p in ps
        if haskey(ic.discrete_idx, p)
            if inner isa Vector{SymbolicT}
                for buf in inner
                end
            end
        end
    end
    if portion isa SciMLStructures.Tunable
    end
end
