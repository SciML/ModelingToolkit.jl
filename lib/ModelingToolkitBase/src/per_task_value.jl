# `PerTaskValue` exists to fix https://github.com/SciML/ModelingToolkit.jl/issues/5336
#
# The per-task mutable state attached to a `System` (the `IRStructure` field and the
# `MutableCacheKey` metadata) used to be stored in `TaskLocalValue`s. A
# `TaskLocalValue` keeps its value in `task_local_storage()` keyed by itself, with a
# strong reference on both the key and the value. Entries are only removed by an
# explicit `empty!` or when the task exits, so every `TaskLocalValue` created — one
# per `System` construction and one per `setproperties` — permanently inserted an
# entry into the current task's storage. Compiling systems in a loop therefore leaked
# ~62 task-local entries per `mtkcompile`, retaining all of their reachable symbols.
#
# `PerTaskValue` inverts the ownership: the per-task values live in a dictionary
# owned by the `System` (or its metadata), keyed weakly on the `Task`. Dropping the
# `System` frees every per-task value, and a dead task's key cannot keep its value
# alive either. `Base.WeakKeyDict` cannot be used for this: its `setindex!` registers
# a finalizer on every key, and the finalizer closure holds the dictionary — i.e.
# the `Task` key would keep the dictionary (and all per-task values) alive for the
# task's lifetime, reproducing the same leak. Hence a plain `Dict{WeakRef, T}` with a
# lock, where entries whose tasks have finished are dropped lazily on write.

"""
    $(TYPEDEF)

A value that is stored once per task which accesses it, initialized lazily by calling
`initializer`. Equivalent to `TaskLocalValue` from TaskLocalValues.jl, except the
per-task values are stored in a `Dict{WeakRef, T}` keyed on the `Task` and owned by
this object rather than in `task_local_storage()`. This means the values are
garbage-collected together with this object (typically a field of `System` or an
entry in its metadata), instead of being retained in the task's local storage until
the task exits.
"""
struct PerTaskValue{T, F}
    initializer::F
    dict::Dict{WeakRef, T}
    lock::Base.Threads.SpinLock
end

function PerTaskValue{T}(initializer::F) where {T, F}
    return PerTaskValue{T, F}(initializer, Dict{WeakRef, T}(), Threads.SpinLock())
end

function _cleanup_dead_tasks!(dict::Dict{WeakRef})
    return filter!(p -> p.first.value !== nothing, dict)
end

Base.eltype(::PerTaskValue{T}) where {T} = T

function Base.getindex(val::PerTaskValue{T}) where {T}
    task = current_task()
    return lock(val.lock) do
        haskey(val.dict, task) && return val.dict[task]
        _cleanup_dead_tasks!(val.dict)
        return val.dict[WeakRef(task)] = val.initializer()::T
    end
end

function Base.setindex!(val::PerTaskValue{T}, value) where {T}
    task = current_task()
    lock(val.lock) do
        _cleanup_dead_tasks!(val.dict)
        val.dict[WeakRef(task)] = convert(T, value)
    end
    return value
end

function Base.empty!(val::PerTaskValue)
    lock(val.lock) do
        delete!(val.dict, current_task())
    end
    return val
end

# The per-task values are lazily created and specific to each `Task` (which cannot be
# serialized or deepcopied). A (de)serialized or deepcopied `PerTaskValue` only retains
# its `initializer` and recreates the per-task values on demand, matching the semantics
# `TaskLocalValue` had.
function Serialization.serialize(s::Serialization.AbstractSerializer, val::PerTaskValue)
    Serialization.writetag(s.io, Serialization.OBJECT_TAG)
    Serialization.serialize(s, typeof(val))
    return Serialization.serialize(s, val.initializer)
end

function Serialization.deserialize(
        s::Serialization.AbstractSerializer, ::Type{PerTaskValue{T, F}}
    ) where {T, F}
    initializer = Serialization.deserialize(s)::F
    return PerTaskValue{T}(initializer)
end

function Base.deepcopy_internal(val::PerTaskValue{T, F}, ::IdDict) where {T, F}
    return PerTaskValue{T}(val.initializer)
end
