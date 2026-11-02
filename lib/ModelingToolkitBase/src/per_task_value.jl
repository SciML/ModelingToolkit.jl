# `PerTaskValue` keeps per-task values in a `Dict{WeakRef, T}` owned by this object
# instead of `task_local_storage()`, so they are released when the owner is collected.
# `Base.WeakKeyDict` is unsuitable for this: its `setindex!` registers a finalizer on
# each key that holds the dictionary, retaining the owner for the task's lifetime.
# Instead, each inserted task gets a `finalizer` holding only a `WeakRef` to this
# object, which sweeps dead entries once the task is collected — so collected tasks
# release their values even while the owner and other tasks keep using it.

"""
    $(TYPEDEF)

A value stored once per `Task` that accesses it, initialized lazily by calling
`initializer`. Equivalent to `TaskLocalValue` from TaskLocalValues.jl, except the
per-task values are stored in a `Dict{WeakRef, T}` owned by this object (typically a
field of `System` or an entry in its metadata) instead of `task_local_storage()`, so
they are released together with this object or once the owning task is collected.
"""
mutable struct PerTaskValue{T, F}
    initializer::F
    dict::Dict{WeakRef, T}
    lock::Base.Threads.SpinLock
end

function PerTaskValue{T}(initializer::F) where {T, F}
    return PerTaskValue{T, F}(initializer, Dict{WeakRef, T}(), Threads.SpinLock())
end

Base.eltype(::PerTaskValue{T}) where {T} = T

function _cleanup_dead_tasks!(dict::Dict{WeakRef})
    return filter!(p -> p.first.value !== nothing, dict)
end

function _sweep_dead_tasks(ref::WeakRef, ::Task)
    val = ref.value
    val === nothing && return nothing
    return lock(val.lock) do
        _cleanup_dead_tasks!(val.dict)
    end
end

# `ref` is the only reference the `finalizer` holds to `val`, so registering it on the
# task cannot keep `val` alive.
function _register_task_cleanup(val::PerTaskValue, task::Task)
    return finalizer(Base.Fix1(_sweep_dead_tasks, WeakRef(val)), task)
end

function Base.getindex(val::PerTaskValue{T}) where {T}
    task = current_task()
    return lock(val.lock) do
        haskey(val.dict, task) && return val.dict[task]
        _cleanup_dead_tasks!(val.dict)
        _register_task_cleanup(val, task)
        return val.dict[WeakRef(task)] = val.initializer()::T
    end
end

function Base.setindex!(val::PerTaskValue{T}, value) where {T}
    task = current_task()
    lock(val.lock) do
        _cleanup_dead_tasks!(val.dict)
        haskey(val.dict, task) || _register_task_cleanup(val, task)
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
# `TaskLocalValue` had. The cycle bookkeeping keeps repeated references to the same
# `PerTaskValue` aliased in the result.
function Serialization.serialize(s::Serialization.AbstractSerializer, val::PerTaskValue)
    Serialization.serialize_cycle_header(s, val) && return nothing
    Serialization.serialize(s, val.initializer)
    return nothing
end

function Serialization.deserialize(
        s::Serialization.AbstractSerializer, ::Type{PerTaskValue{T, F}}
    ) where {T, F}
    initializer = Serialization.deserialize(s)::F
    val = PerTaskValue{T}(initializer)
    Serialization.deserialize_cycle(s, val)
    return val
end

function Base.deepcopy_internal(val::PerTaskValue{T, F}, iddict::IdDict) where {T, F}
    haskey(iddict, val) && return iddict[val]::PerTaskValue{T, F}
    return iddict[val] = PerTaskValue{T}(val.initializer)
end
