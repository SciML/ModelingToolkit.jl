"""
    $TYPEDEF

A struct which stores dependency information about the bound parameters in a system.

# Fields

$TYPEDFIELDS
"""
struct ParameterBindingsGraph
    """
    An ordered set of bound parameters, in topologically sorted order.
    """
    bound_ps::AtomicArraySet{OrderedDict{SymbolicT, Nothing}}
    """
    Since indexing `OrderedSet` is deprecated, this maps an index to the corresponding
    bound parameter. Literally just `collect(bound_ps)`.
    """
    bound_ps_order::Vector{SymbolicT}
    """
    Mapping from bound parameters to their index in `bound_ps`.
    """
    bound_par_idx::AtomicArrayDict{Int, Dict{SymbolicT, Int}}
    """
    Dependency graph for bindings. Edges flow from a bound parameter to the
    other bound parameters it depends on.
    """
    dependency_graph::SimpleDiGraph{Int}
end

function ParameterBindingsGraph(sys::AbstractSystem)
    if !isempty(get_systems(sys))
        throw(ArgumentError("`ParameterBindingsGraph` can only be created from a flattened system."))
    end
    all_ps = AtomicArraySet{OrderedDict{SymbolicT, Nothing}}()
    for p in get_ps(sys)
        push_as_atomic_array!(all_ps, p)
    end
    # If the system already has a pbgraph, it was previously `complete`d. Those bound parameters
    # will not be in `ps`, but need to be considered here.
    old_pbgraph = get_parameter_bindings_graph(sys)
    if old_pbgraph isa ParameterBindingsGraph
        union!(all_ps, old_pbgraph.bound_ps)
    end
    # This may be called on a non-completed system, or a `split=false` system, or one with
    # callbacks modified. The only fool-proof way to find discretes is to filter out ones
    # from callbacks.
    all_discretes = get_all_discretes(sys)

    static_ps = setdiff!(all_ps, all_discretes)
    binds = bindings(sys)
    bound_ps = intersect!(static_ps, keys(binds))
    filter!(!(Base.Fix2(===, COMMON_MISSING) ∘ Base.Fix1(getindex, binds)), bound_ps)

    bound_par_idx = AtomicArrayDict{Int}()
    for (i, p) in enumerate(bound_ps)
        bound_par_idx[p] = i
    end

    dependency_graph = SimpleDiGraph{Int}(length(bound_ps))
    varsbuf = Set{SymbolicT}()
    for (i, p) in enumerate(bound_ps)
        val = binds[p]
        empty!(varsbuf)
        Symbolics.get_variables!(varsbuf, val, bound_ps)

        for dep in varsbuf
            add_edge!(dependency_graph, i, bound_par_idx[dep])
        end
    end

    possible_cycles = simplecycles_iter(dependency_graph, 1)
    if !isempty(possible_cycles)
        throw(CyclicBindingsError(collect(bound_ps)[possible_cycles[1]]))
    end

    toporder = topological_sort(dependency_graph)
    reverse!(toporder)
    bound_ps_order = collect(bound_ps)
    _bound_ps = AtomicArraySet{OrderedDict{SymbolicT, Nothing}}()
    _bound_par_idx = AtomicArrayDict{Int}()
    _dep_graph = SimpleDiGraph{Int}(length(bound_ps))

    for (ii, i) in enumerate(toporder)
        push!(_bound_ps, bound_ps_order[i])
        _bound_par_idx[bound_ps_order[i]] = ii
    end

    for e in edges(dependency_graph)
        add_edge!(_dep_graph, _bound_par_idx[bound_ps_order[src(e)]], _bound_par_idx[bound_ps_order[dst(e)]])
    end
    bound_ps = _bound_ps
    bound_par_idx = _bound_par_idx
    dependency_graph = _dep_graph

    return ParameterBindingsGraph(bound_ps, collect(bound_ps), bound_par_idx, dependency_graph)
end

"""
    $TYPEDSIGNATURES

Find the bound parameters used by `expr`.
"""
function bound_parameters_used_by!(buffer::OrderedSet{SymbolicT}, sys::AbstractSystem, expr;
                                   bgraph::ParameterBindingsGraph = get_parameter_bindings_graph(sys))
    # No point searching if we've already included all of them.
    if issubset(bgraph.bound_ps, buffer)
        return buffer
    end

    vars = OrderedSet{SymbolicT}()
    Symbolics.get_variables!(vars, expr, bgraph.bound_ps)
    idxs = Int[]
    for v in vars
        push!(idxs, bgraph.bound_par_idx[v])
    end
    for i in BFSIterator(bgraph.dependency_graph, idxs)
        push!(buffer, bgraph.bound_ps_order[i])
    end

    return buffer
end

"""
    $TYPEDSIGNATURES

Topologically sort the bound parameters `bound_ps`.
"""
function sort_bound_parameters!(bound_ps::OrderedSet{SymbolicT}, sys::AbstractSystem;
                                bgraph::ParameterBindingsGraph = get_parameter_bindings_graph(sys))
    sort!(bound_ps; by = Base.Fix1(getindex, bgraph.bound_par_idx))
    return bound_ps
end

struct CyclicBindingsError <: Exception
    cycle_ps::Vector{SymbolicT}
end

function Base.showerror(io::IO, err::CyclicBindingsError)
    println(io, """
    The bindings for parameters were found to have at least one cycle involving the \
    follow parameters:
    """)
    for p in err.cycle_ps
        println(io, p)
    end
end
