struct ParameterIndex{P, I}
end
struct IndexCache
end
function IndexCache(sys::AbstractSystem)
    return IndexCache()
end
