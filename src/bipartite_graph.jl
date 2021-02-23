module BipartiteGraphs

export BipartiteEdge, BipartiteGraph

export 𝑠vertices, 𝑑vertices, has_𝑠vertex, has_𝑑vertex, 𝑠neighbors, 𝑑neighbors,
       𝑠edges, 𝑑edges, nsrcs, ndsts, SRC, DST

using DocStringExtensions
using Reexport
using UnPack
using SparseArrays
@reexport using LightGraphs
using Setfield

###
### Edges & Vertex
###
@enum VertType SRC DST ALL

struct BipartiteEdge{I<:Integer} <: LightGraphs.AbstractEdge{I}
    src::I
    dst::I
    function BipartiteEdge(src::I, dst::V) where {I,V}
        T = promote_type(I, V)
        new{T}(T(src), T(dst))
    end
end

LightGraphs.src(edge::BipartiteEdge) = edge.src
LightGraphs.dst(edge::BipartiteEdge) = edge.dst

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
mutable struct BipartiteGraph{I<:Integer,M} <: LightGraphs.AbstractGraph{I}
    ne::Int
    fadjlist::Vector{Vector{I}} # `fadjlist[src] => dsts`
    badjlist::Vector{Vector{I}} # `badjlist[dst] => srcs`
    metadata::M
end
BipartiteGraph(ne::Integer, fadj::AbstractVector, badj::AbstractVector) = BipartiteGraph(ne, fadj, badj, nothing)

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
function BipartiteGraph(nsrcs::T, ndsts::T; metadata=nothing) where T
    fadjlist = map(_->T[], 1:nsrcs)
    badjlist = map(_->T[], 1:ndsts)
    BipartiteGraph(0, fadjlist, badjlist, metadata)
end

Base.eltype(::Type{<:BipartiteGraph{I}}) where I = I
function Base.empty!(g::BipartiteGraph)
    foreach(empty!, g.fadjlist)
    foreach(empty!, g.badjlist)
    g.ne = 0
    if g.metadata !== nothing
        foreach(empty!, g.metadata)
    end
    g
end
Base.length(::BipartiteGraph) = error("length is not well defined! Use `ne` or `nv`.")

if isdefined(LightGraphs, :has_contiguous_vertices)
    LightGraphs.has_contiguous_vertices(::Type{<:BipartiteGraph}) = false
end
LightGraphs.is_directed(::Type{<:BipartiteGraph}) = false
LightGraphs.vertices(g::BipartiteGraph) = (𝑠vertices(g), 𝑑vertices(g))
𝑠vertices(g::BipartiteGraph) = axes(g.fadjlist, 1)
𝑑vertices(g::BipartiteGraph) = axes(g.badjlist, 1)
has_𝑠vertex(g::BipartiteGraph, v::Integer) = v in 𝑠vertices(g)
has_𝑑vertex(g::BipartiteGraph, v::Integer) = v in 𝑑vertices(g)
𝑠neighbors(g::BipartiteGraph, i::Integer, with_metadata::Val{M}=Val(false)) where M = M ? zip(g.fadjlist[i], g.metadata[i]) : g.fadjlist[i]
𝑑neighbors(g::BipartiteGraph, j::Integer, with_metadata::Val{M}=Val(false)) where M = M ? zip(g.badjlist[j], (g.metadata[i][j] for i in g.badjlist[j])) : g.badjlist[j]
LightGraphs.ne(g::BipartiteGraph) = g.ne
LightGraphs.nv(g::BipartiteGraph) = sum(length, vertices(g))
LightGraphs.edgetype(g::BipartiteGraph{I}) where I = BipartiteEdge{I}

nsrcs(g::BipartiteGraph) = length(𝑠vertices(g))
ndsts(g::BipartiteGraph) = length(𝑑vertices(g))

function LightGraphs.has_edge(g::BipartiteGraph, edge::BipartiteEdge)
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

LightGraphs.add_edge!(g::BipartiteGraph, i::Integer, j::Integer, md=NO_METADATA) = add_edge!(g, BipartiteEdge(i, j), md)
function LightGraphs.add_edge!(g::BipartiteGraph, edge::BipartiteEdge, md=NO_METADATA)
    @unpack fadjlist, badjlist = g
    verts = vertices(g)
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
    @inbounds list = badjlist[d]
    index = searchsortedfirst(list, s)
    insert!(list, index, s)
    return true  # edge successfully added
end

function LightGraphs.add_vertex!(g::BipartiteGraph{T}, type::VertType) where T
    if type === DST
        push!(g.badjlist, T[])
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
LightGraphs.edges(g::BipartiteGraph) = BipartiteEdgeIter(g, Val(ALL))
𝑠edges(g::BipartiteGraph) = BipartiteEdgeIter(g, Val(SRC))
𝑑edges(g::BipartiteGraph) = BipartiteEdgeIter(g, Val(DST))

struct BipartiteEdgeIter{T,G} <: LightGraphs.AbstractEdgeIter
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
function LightGraphs.incidence_matrix(g::BipartiteGraph, val=true)
    I = Int[]
    J = Int[]
    for i in 𝑠vertices(g), n in 𝑠neighbors(g, i)
        push!(I, i)
        push!(J, n)
    end
    S = sparse(I, J, val, nsrcs(g), ndsts(g))
end

end # module
