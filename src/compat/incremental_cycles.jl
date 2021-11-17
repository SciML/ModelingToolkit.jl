using Base.Iterators: repeated

# Abstract Interface

"""
    abstract type IncrementalCycleTracker

The supertype for incremental cycle detection problems. The abstract type
constructor IncrementalCycleTracker(G) may be used to automatically select
a specific incremental cycle detection algorithm. See [`add_edge_checked!`](@ref)
for a usage example.
"""
abstract type IncrementalCycleTracker{I} <: AbstractGraph{I} end

function (::Type{IncrementalCycleTracker})(s::AbstractGraph{I}; dir=:out) where {I}
    # TODO: Once we have more algorithms, the poly-algorithm decision goes here.
    # For now, we only have Algorithm N.
    return DenseGraphICT_BFGT_N{something(dir == :in, false)}(s)
end

# Cycle Detection Interface
"""
    add_edge_checked!([f!,], ict::IncrementalCycleTracker, v, w)

Using the incremental cycle tracker, ict, check whether adding the edge `v=>w`.
Would introduce a cycle in the underlying graph. If so, return false and leave
the ict intact. If not, update the underlying graph and return true.

# Optional `f!` Argument

By default the `add_edge!` function is used to update the underlying graph.
However, for more complicated graphs, users may wish to manually specify the
graph update operation. This may be accomplished by passing the optional `f!`
callback arhgument. This callback is called on the underlying graph when no
cycle is detected and is required to modify the underlying graph in order to
effectuate the proposed edge addition.

# Batched edge additions

Optionally, either `v` or `w` (depending on the `in_out_reverse` flag) may be a
collection of vertices representing a batched addition of vertices sharing a
common source or target more efficiently than individual updates.

## Example

```jldoctest
julia> G = SimpleDiGraph(3)

julia> ict = IncrementalCycleTracker(G)
BFGT_N cycle tracker on {3, 0} directed simple Int64 graph

julia> add_edge_checked!(ict, 1, 2)
true

julia> collect(edges(G))
1-element Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}:
 Edge 1 => 2

julia> add_edge_checked!(ict, 2, 3)
true

julia> collect(edges(G))
2-element Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}:
 Edge 1 => 2
 Edge 2 => 3

julia> add_edge_checked!(ict, 3, 1) # Would add a cycle
false

julia> collect(edges(G))
2-element Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}:
Edge 1 => 2
Edge 2 => 3
```
"""
function add_edge_checked! end

to_edges(v::Integer, w::Integer) = (v=>w,)
to_edges(v::Integer, ws) = zip(repeated(v), ws)
to_edges(vs, w::Integer) = zip(vs, repeated(w))

add_edge_checked!(ict::IncrementalCycleTracker, vs, ws) = add_edge_checked!(ict, vs, ws) do g
    foreach(((v, w),)->add_edge!(g, v, w), to_edges(vs, ws))
end

# Utilities
"""
    struct TransactionalVector

A vector with one checkpoint that may be reverted to by calling `revert!`. The setpoint itself
is set by calling `commit!`.
"""
struct TransactionalVector{T} <: AbstractVector{T}
    v::Vector{T}
    log::Vector{Pair{Int, T}}
    TransactionalVector(v::Vector{T}) where {T} =
        new{T}(v, Vector{Pair{Int, T}}())
end

function commit!(v::TransactionalVector)
    empty!(v.log)
    return nothing
end

function revert!(vec::TransactionalVector)
    for (idx, val) in reverse(vec.log)
        vec.v[idx] = val
    end
    return nothing
end

function Base.setindex!(vec::TransactionalVector, val, idx)
    oldval = vec.v[idx]
    vec.v[idx] = val
    push!(vec.log, idx=>oldval)
    return nothing
end
Base.getindex(vec::TransactionalVector, idx) = vec.v[idx]
Base.size(vec::TransactionalVector) = size(vec.v)

# Specific Algorithms

const bibliography = """
## References

[BFGT15] Michael A. Bender, Jeremy T. Fineman, Seth Gilbert, and Robert E. Tarjan. 2015
    A New Approach to Incremental Cycle Detection and Related Problems.
    ACM Trans. Algorithms 12, 2, Article 14 (December 2015), 22 pages.
    DOI: http://dx.doi.org/10.1145/2756553
"""

## Bender, Algorithm N

"""
    struct DenseGraphICT_BFGT_N

Implements the "Naive" (Algorithm N) Bender-Fineman-Gilbert-Tarjan one-way line search incremental cycle detector
for dense graphs from [BFGT15] (Section 3).

$bibliography
"""
struct DenseGraphICT_BFGT_N{InOutReverse, I, G<:AbstractGraph{I}} <: IncrementalCycleTracker{I}
    graph::G
    levels::TransactionalVector{Int}
    DenseGraphICT_BFGT_N{InOutReverse}(g::G) where {InOutReverse, I, G<:AbstractGraph{I}} =
        new{InOutReverse, I, G}(g, TransactionalVector(fill(0, nv(g))))
end
function Base.show(io::IO, ict::DenseGraphICT_BFGT_N)
    print(io, "BFGT_N cycle tracker on ")
    show(io, ict.graph)
end

function topological_sort(ict::DenseGraphICT_BFGT_N{InOutReverse}) where {InOutReverse}
    # The ICT levels are a weak topological ordering, so a sort of the levels
    # will give a topological sort of the vertices.
    perm = sortperm(ict.levels)
    InOutReverse && (perm = reverse(perm))
    return perm
end

# Even when both `v` and `w` are integer, we know that `v` would come first, so
# we prefer to check for `v` as the cycle vertex in this case.
add_edge_checked!(f!, ict::DenseGraphICT_BFGT_N{false}, v::Integer, ws) =
    _check_cycle_add!(f!, ict, to_edges(v, ws), v)
add_edge_checked!(f!, ict::DenseGraphICT_BFGT_N{true}, vs, w::Integer) =
    _check_cycle_add!(f!, ict, to_edges(vs, w), w)

### [BFGT15] Algorithm N
#
# Implementation Notes
#
# This is Algorithm N from [BFGT15] (Section 3), plus limited patching support and
# a number of standard tricks. Namely:
#
# 1. Batching is supported as long as there is only a single source or destination
#    vertex. General batching is left as an open problem. The reason that the
#    single source/dest batching is easy to add is that we know that either the
#    source or the destination vertex is guaranteed to be a part of any cycle
#    that we may have added. Thus we're guaranteed to encounter one of the two
#    verticies in our cycle validation and the rest of the algorithm goes through
#    as usual.
# 2. We opportunistically traverse each edge when we see it and only add it
#    to the worklist if we know that traversal will recurse further.
# 3. We add some early out checks to detect we're about to do redundant work.
function _check_cycle_add!(f!, ict::DenseGraphICT_BFGT_N{InOutReverse}, edges, v) where {InOutReverse}
    g = ict.graph
    worklist = Pair{Int, Int}[]
    # TODO: In the case where there's a single target vertex, we could saturate
    # the level first before we assign it to the tracked vector to save some
    # log space.
    for (v, w) in edges
        InOutReverse && ((v, w) = (w, v))
        if ict.levels[v] < ict.levels[w]
            continue
        end
        v == w && return false
        ict.levels[w] = ict.levels[v] + 1
        push!(worklist, v=>w)
    end
    while !isempty(worklist)
        (x, y) = popfirst!(worklist)
        xlevel = ict.levels[x]
        ylevel = ict.levels[y]
        if xlevel >= ylevel
            # The xlevel may have been incremented further since we added this
            # edge to the worklist.
            ict.levels[y] = ylevel = xlevel + 1
        elseif ylevel > xlevel + 1
            # Some edge traversal scheduled for later already incremented this
            # level past where we would have been. Delay processing until then.
            continue
        end
        for z in (InOutReverse ? inneighbors(g, y) : outneighbors(g, y))
            if z == v
                revert!(ict.levels)
                return false
            end
            if ylevel >= ict.levels[z]
                ict.levels[z] = ylevel + 1
                push!(worklist, y=>z)
            end
        end
    end
    commit!(ict.levels)
    f!(g)
    return true
end
