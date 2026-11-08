# `PerTaskValue` keeps per-task values in a `WeakKeyDict{Task, T}` owned by this object
# instead of `task_local_storage()`: entries are keyed weakly on the `Task`, so values
# are released once the owning task is collected or this object is.

"""
    $(TYPEDEF)

A value stored once per `Task` that accesses it, initialized lazily by calling
`initializer`. Equivalent to `TaskLocalValue` from TaskLocalValues.jl, except the
per-task values are stored in a `WeakKeyDict{Task, T}` owned by this object (typically
a field of `System` or an entry in its metadata) instead of `task_local_storage()`, so
they are released together with this object or once the owning task is collected.
"""
struct PerTaskValue{T, F}
    initializer::F
    values::WeakKeyDict{Task, T}
end

function PerTaskValue{T}(initializer::F) where {T, F}
    return PerTaskValue{T, F}(initializer, WeakKeyDict{Task, T}())
end

Base.eltype(::PerTaskValue{T}) where {T} = T

function Base.getindex(val::PerTaskValue{T}) where {T}
    return get!(val.initializer, val.values, current_task())::T
end

function Base.setindex!(val::PerTaskValue{T}, value) where {T}
    val.values[current_task()] = convert(T, value)
    return value
end

function Base.empty!(val::PerTaskValue)
    delete!(val.values, current_task())
    return val
end
