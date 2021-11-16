module BipartiteGraphs

export BipartiteEdge, BipartiteGraph, DiCMOBiGraph, Unassigned, unassigned

export 𝑠vertices, 𝑑vertices, has_𝑠vertex, has_𝑑vertex, 𝑠neighbors, 𝑑neighbors,
       𝑠edges, 𝑑edges, nsrcs, ndsts, SRC, DST

using DocStringExtensions
using UnPack
using SparseArrays
using Graphs
using Setfield

### Matching
struct Unassigned
    global unassigned
    const unassigned = Unassigned.instance
end

###
### Edges & Vertex
###
@enum VertType SRC DST ALL

struct BipartiteEdge{I<:Integer} <: Graphs.AbstractEdge{I}
    src::I
    dst::I
    function BipartiteEdge(src::I, dst::V) where {I,V}
        T = promote_type(I, V)
        new{T}(T(src), T(dst))
    end
end

Graphs.src(edge::BipartiteEdge) = edge.src
Graphs.dst(edge::BipartiteEdge) = edge.dst

function Base.show(io::IO, edge::BipartiteEdge)
    @unpack src, dst = edge
    print(io, "[src: ", src, "] => [dst: ", dst, "]")
end

Base.:(==)(a::BipartiteEdge, b::BipartiteEdge) = src(a) == src(b) && dst(a) == dst(b)

###
### Graph
###
"""
$(TYPEDEF)

A bipartite graph representation between two, possibly distinct, sets of vertices
(source and dependencies). Maps source vertices, labelled `1:N₁`, to vertices
on which they depend (labelled `1:N₂`).

# Fields
$(FIELDS)

# Example
```julia
using ModelingToolkit

ne = 4
srcverts = 1:4
depverts = 1:2

# six source vertices
fadjlist = [[1],[1],[2],[2],[1],[1,2]]

# two vertices they depend on
badjlist = [[1,2,5,6],[3,4,6]]

bg = BipartiteGraph(7, fadjlist, badjlist)
```
"""
mutable struct BipartiteGraph{I<:Integer, M} <: LightGraphs.AbstractGraph{I}
    ne::Int
    fadjlist::Vector{Vector{I}} # `fadjlist[src] => dsts`
    badjlist::Union{Vector{Vector{I}},I} # `badjlist[dst] => srcs` or `ndsts`
    metadata::M
end
BipartiteGraph(ne::Integer, fadj::AbstractVector, badj::Union{AbstractVector,Integer}=maximum(maximum, fadj); metadata=nothing) = BipartiteGraph(ne, fadj, badj, metadata)

"""
```julia
Base.isequal(bg1::BipartiteGraph{T}, bg2::BipartiteGraph{T}) where {T<:Integer}
```

Test whether two [`BipartiteGraph`](@ref)s are equal.
"""
function Base.isequal(bg1::BipartiteGraph{T}, bg2::BipartiteGraph{T}) where {T<:Integer}
    iseq = (bg1.ne == bg2.ne)
    iseq &= (bg1.fadjlist == bg2.fadjlist)
    iseq &= (bg1.badjlist == bg2.badjlist)
    iseq
end

"""
$(SIGNATURES)

Build an empty `BipartiteGraph` with `nsrcs` sources and `ndsts` destinations.
"""
function BipartiteGraph(nsrcs::T, ndsts::T, backedge::Val{B}=Val(true); metadata=nothing) where {T,B}
    fadjlist = map(_->T[], 1:nsrcs)
    badjlist = B ? map(_->T[], 1:ndsts) : ndsts
    BipartiteGraph(0, fadjlist, badjlist, metadata)
end

Base.eltype(::Type{<:BipartiteGraph{I}}) where I = I
function Base.empty!(g::BipartiteGraph)
    foreach(empty!, g.fadjlist)
    g.badjlist isa AbstractVector && foreach(empty!, g.badjlist)
    g.ne = 0
    if g.metadata !== nothing
        foreach(empty!, g.metadata)
    end
    g
end
Base.length(::BipartiteGraph) = error("length is not well defined! Use `ne` or `nv`.")

@noinline throw_no_back_edges() = throw(ArgumentError("The graph has no back edges."))

if isdefined(Graphs, :has_contiguous_vertices)
    Graphs.has_contiguous_vertices(::Type{<:BipartiteGraph}) = false
end
Graphs.is_directed(::Type{<:BipartiteGraph}) = false
Graphs.vertices(g::BipartiteGraph) = (𝑠vertices(g), 𝑑vertices(g))
𝑠vertices(g::BipartiteGraph) = axes(g.fadjlist, 1)
𝑑vertices(g::BipartiteGraph) = g.badjlist isa AbstractVector ? axes(g.badjlist, 1) : Base.OneTo(g.badjlist)
has_𝑠vertex(g::BipartiteGraph, v::Integer) = v in 𝑠vertices(g)
has_𝑑vertex(g::BipartiteGraph, v::Integer) = v in 𝑑vertices(g)
𝑠neighbors(g::BipartiteGraph, i::Integer, with_metadata::Val{M}=Val(false)) where M = M ? zip(g.fadjlist[i], g.metadata[i]) : g.fadjlist[i]
function 𝑑neighbors(g::BipartiteGraph, j::Integer, with_metadata::Val{M}=Val(false)) where M
    g.badjlist isa AbstractVector || throw_no_back_edges()
    M ? zip(g.badjlist[j], (g.metadata[i][j] for i in g.badjlist[j])) : g.badjlist[j]
end
Graphs.ne(g::BipartiteGraph) = g.ne
Graphs.nv(g::BipartiteGraph) = sum(length, vertices(g))
Graphs.edgetype(g::BipartiteGraph{I}) where I = BipartiteEdge{I}

nsrcs(g::BipartiteGraph) = length(𝑠vertices(g))
ndsts(g::BipartiteGraph) = length(𝑑vertices(g))

function Graphs.has_edge(g::BipartiteGraph, edge::BipartiteEdge)
    @unpack src, dst = edge
    (src in 𝑠vertices(g) && dst in 𝑑vertices(g)) || return false  # edge out of bounds
    insorted(𝑠neighbors(src), dst)
end

###
### Populate
###
struct NoMetadata
end
const NO_METADATA = NoMetadata()

Graphs.add_edge!(g::BipartiteGraph, i::Integer, j::Integer, md=NO_METADATA) = add_edge!(g, BipartiteEdge(i, j), md)
function Graphs.add_edge!(g::BipartiteGraph, edge::BipartiteEdge, md=NO_METADATA)
    @unpack fadjlist, badjlist = g
    s, d = src(edge), dst(edge)
    (has_𝑠vertex(g, s) && has_𝑑vertex(g, d)) || error("edge ($edge) out of range.")
    @inbounds list = fadjlist[s]
    index = searchsortedfirst(list, d)
    @inbounds (index <= length(list) && list[index] == d) && return false  # edge already in graph
    insert!(list, index, d)
    if md !== NO_METADATA
        insert!(g.metadata[s], index, md)
    end

    g.ne += 1
    if badjlist isa AbstractVector
        @inbounds list = badjlist[d]
        index = searchsortedfirst(list, s)
        insert!(list, index, s)
    end
    return true  # edge successfully added
end

function Graphs.add_vertex!(g::BipartiteGraph{T}, type::VertType) where T
    if type === DST
        if g.badjlist isa AbstractVector
            push!(g.badjlist, T[])
        else
            g.badjlist += 1
        end
    elseif type === SRC
        push!(g.fadjlist, T[])
    else
        error("type ($type) must be either `DST` or `SRC`")
    end
    return true  # vertex successfully added
end

###
### Edges iteration
###
Graphs.edges(g::BipartiteGraph) = BipartiteEdgeIter(g, Val(ALL))
𝑠edges(g::BipartiteGraph) = BipartiteEdgeIter(g, Val(SRC))
𝑑edges(g::BipartiteGraph) = BipartiteEdgeIter(g, Val(DST))

struct BipartiteEdgeIter{T,G} <: Graphs.AbstractEdgeIter
    g::G
    type::Val{T}
end

Base.length(it::BipartiteEdgeIter) = ne(it.g)
Base.length(it::BipartiteEdgeIter{ALL}) = 2ne(it.g)

Base.eltype(it::BipartiteEdgeIter) = edgetype(it.g)

function Base.iterate(it::BipartiteEdgeIter{SRC,<:BipartiteGraph{T}}, state=(1, 1, SRC)) where T
    @unpack g = it
    neqs = nsrcs(g)
    neqs == 0 && return nothing
    eq, jvar = state

    while eq <= neqs
        eq′ = eq
        vars = 𝑠neighbors(g, eq′)
        if jvar > length(vars)
            eq += 1
            jvar = 1
            continue
        end
        edge = BipartiteEdge(eq′, vars[jvar])
        state = (eq, jvar + 1, SRC)
        return edge, state
    end
    return nothing
end

function Base.iterate(it::BipartiteEdgeIter{DST,<:BipartiteGraph{T}}, state=(1, 1, DST)) where T
    @unpack g = it
    nvars = ndsts(g)
    nvars == 0 && return nothing
    ieq, jvar = state

    while jvar <= nvars
        eqs = 𝑑neighbors(g, jvar)
        if ieq > length(eqs)
            ieq = 1
            jvar += 1
            continue
        end
        edge = BipartiteEdge(eqs[ieq], jvar)
        state = (ieq + 1, jvar, DST)
        return edge, state
    end
    return nothing
end

function Base.iterate(it::BipartiteEdgeIter{ALL,<:BipartiteGraph}, state=nothing)
    if state === nothing
        ss = iterate((@set it.type = Val(SRC)))
    elseif state[3] === SRC
        ss = iterate((@set it.type = Val(SRC)), state)
    elseif state[3] == DST
        ss = iterate((@set it.type = Val(DST)), state)
    end
    if ss === nothing && state[3] == SRC
        return iterate((@set it.type = Val(DST)))
    else
        return ss
    end
end

###
### Utils
###
function Graphs.incidence_matrix(g::BipartiteGraph, val=true)
    I = Int[]
    J = Int[]
    for i in 𝑠vertices(g), n in 𝑠neighbors(g, i)
        push!(I, i)
        push!(J, n)
    end
    S = sparse(I, J, val, nsrcs(g), ndsts(g))
end


"""
    struct DiCMOBiGraph

This data structure implements a "directed, contracted, matching-oriented" view of an
original (undirected) bipartite graph. It has two modes, depending on the `Transposed`
flag, which switches the direction of the induced matching.

Essentially the graph adapter performs two largely orthogonal functions
[`Transposed == true` differences are indicated in square brackets]:

1. It pairs an undirected bipartite graph with a matching of the destination vertex.

    This matching is used to induce an orientation on the otherwise undirected graph:
    Matched edges pass from destination to source [source to desination], all other edges
    pass in the opposite direction.

2. It exposes the graph view obtained by contracting the destination [source] vertices
   along the matched edges.

The result of this operation is an induced, directed graph on the source [destination] vertices.
The resulting graph has a few desirable properties. In particular, this graph
is acyclic if and only if the induced directed graph on the original bipartite
graph is acyclic.

"""
struct DiCMOBiGraph{Transposed, I, G<:BipartiteGraph{I}, M} <: Graphs.AbstractGraph{I}
    graph::G
    matching::M
    DiCMOBiGraph{Transposed}(g::G, m::M) where {Transposed, I, G<:BipartiteGraph{I}, M} =
        new{Transposed, I, G, M}(g, m)
end
DiCMOBiGraph{Transposed}(g::BipartiteGraph) where {Transposed} = DiCMOBiGraph{Transposed}(g, Union{Unassigned, Int}[unassigned for i = 1:ndsts(g)])
Graphs.is_directed(::Type{<:DiCMOBiGraph}) = true
Graphs.nv(g::DiCMOBiGraph{Transposed}) where {Transposed} = Transposed ? ndsts(g.graph) : nsrcs(g.graph)
Graphs.vertices(g::DiCMOBiGraph{Transposed}) where {Transposed} = Transposed ? 𝑑vertices(g.graph) : 𝑠vertices(g.graph)

struct CMONeighbors{Transposed, V}
    g::DiCMOBiGraph{Transposed}
    v::V
    CMONeighbors{Transposed}(g::DiCMOBiGraph{Transposed}, v::V) where {Transposed, V} =
        new{Transposed, V}(g, v)
end
Graphs.outneighbors(g::DiCMOBiGraph{false}, v) = CMONeighbors{false}(g, v)
Base.iterate(c::CMONeighbors{false}) = iterate(c, (c.g.graph.fadjlist[c.v],))
function Base.iterate(c::CMONeighbors{false}, (l, state...))
    while true
        r = iterate(l, state...)
        r === nothing && return nothing
        # If this is a matched edge, skip it, it's reversed in the induced
        # directed graph. Otherwise, if there is no matching for this destination
        # edge, also skip it, since it got delted in the contraction.
        vsrc = c.g.matching[r[1]]
        if vsrc === c.v || vsrc === unassigned
            state = (r[2],)
            continue
        end
        return vsrc, (l, r[2])
    end
end

Graphs.inneighbors(g::DiCMOBiGraph{true}, v) = CMONeighbors{true}(g, v)
function Base.iterate(c::CMONeighbors{true})
    vsrc = c.g.matching[c.v]
    vsrc === unassigned && return nothing
    iterate(c, (c.g.graph.fadjlist[vsrc],))
end
function Base.iterate(c::CMONeighbors{true}, (l, state...))
    while true
        r = iterate(l, state...)
        r === nothing && return nothing
        if r[1] === c.v
            state = (r[2],)
            continue
        end
        return r[1], (l, r[2])
    end
end

end # module
