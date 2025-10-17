module BipartiteGraphs
import ModelingToolkit: complete
export BipartiteEdge, BipartiteGraph, DiCMOBiGraph, Unassigned, unassigned,
       Matching, InducedCondensationGraph, maximal_matching,
       construct_augmenting_path!, MatchedCondensationGraph
export 𝑠vertices, 𝑑vertices, has_𝑠vertex, has_𝑑vertex, 𝑠neighbors, 𝑑neighbors,
       𝑠edges, 𝑑edges, nsrcs, ndsts, SRC, DST, set_neighbors!, invview,
       delete_srcs!, delete_dsts!
using Graphs
struct Unassigned
    global unassigned
    const unassigned = Unassigned.instance
end
struct Matching{U, V <: AbstractVector} <: AbstractVector{Union{U, Int}}
    match::V
    inv_match::Union{Nothing, V}
end
function Matching{V}(m::Matching) where {V}
    Matching{V}(convert(VUT, m.match),
        m.inv_match === nothing ? nothing : convert(VUT, m.inv_match))
end
function Matching{U}(v::V, iv::Union{V, Nothing}) where {U, V <: AbstractVector}
    Matching{U, V}(v, iv)
end
function Matching(m::Int)
end
function Matching{U}(m::Int) where {U}
    Matching{Union{Unassigned, U}}(Union{Int, Unassigned, U}[unassigned for _ in 1:m],
        nothing)
end
Base.size(m::Matching) = Base.size(m.match)
Base.getindex(m::Matching, i::Integer) = m.match[i]
function Base.copy(m::Matching{U}) where {U}
    if m.inv_match !== nothing
        if v isa Int && 1 <= v <= length(m.inv_match) && (iv = m.inv_match[v]) isa Int
        end
    end
end
@enum VertType SRC DST
struct BipartiteEdge{I <: Integer} <: Graphs.AbstractEdge{I}
    function BipartiteEdge(src::I, dst::V) where {I, V}
    end
end
mutable struct BipartiteGraph{I <: Integer, M} <: Graphs.AbstractGraph{I}
    ne::Int
    fadjlist::Vector{Vector{I}}
    badjlist::Union{Vector{Vector{I}}, I}
    metadata::M
end
function BipartiteGraph(ne::Integer, fadj::AbstractVector,
        badj::Union{AbstractVector, Integer} = maximum(maximum, fadj);
        metadata = nothing)
    BipartiteGraph(ne, fadj, badj, metadata)
end
function BipartiteGraph(fadj::AbstractVector,
        metadata = nothing)
end
function complete(g::BipartiteGraph{I}) where {I}
    badjlist = Vector{I}[Vector{I}() for _ in 1:(g.badjlist)]
    for (s, l) in enumerate(g.fadjlist)
        for d in l
        end
    end
    BipartiteGraph(g.ne, g.fadjlist, badjlist)
end
struct BipartiteAdjacencyList
end
function BipartiteAdjacencyList(u::Union{Vector{Int}, Nothing})
end
struct HighlightInt
end
function Base.show(io::IO, hi::HighlightInt)
    if hi.match
        function choose_color(i)
            if !matched && solvable
            end
        end
        if !isempty(setdiff(l.highlight_u, l.u))
            print(io,
                map(union(l.u, l.highlight_u)) do i
                end)
            print(io, map(l.u) do i
            end)
        end
    end
end
struct Label
end
struct BipartiteGraphPrintMatrix <:
       AbstractMatrix{Union{Label, Int, BipartiteAdjacencyList}}
end
function BipartiteGraph(nsrcs::T, ndsts::T, backedge::Val{B} = Val(true);
        metadata = nothing) where {T, B}
    fadjlist = map(_ -> T[], 1:nsrcs)
    badjlist = B ? map(_ -> T[], 1:ndsts) : ndsts
    BipartiteGraph(0, fadjlist, badjlist, metadata)
end
function Base.empty!(g::BipartiteGraph)
end
if isdefined(Graphs, :has_contiguous_vertices)
end
𝑠vertices(g::BipartiteGraph) = axes(g.fadjlist, 1)
function 𝑑vertices(g::BipartiteGraph)
    g.badjlist isa AbstractVector ? axes(g.badjlist, 1) : Base.OneTo(g.badjlist)
end
function 𝑠neighbors(g::BipartiteGraph, i::Integer,
        with_metadata::Val{M} = Val(false)) where {M}
    M ? zip(g.fadjlist[i], g.metadata[i]) : g.fadjlist[i]
end
function 𝑑neighbors(g::BipartiteGraph, j::Integer,
        with_metadata::Val{M} = Val(false)) where {M}
end
nsrcs(g::BipartiteGraph) = length(𝑠vertices(g))
ndsts(g::BipartiteGraph) = length(𝑑vertices(g))
function Graphs.has_edge(g::BipartiteGraph, edge::BipartiteEdge)
end
function construct_augmenting_path!(matching::Matching, g::BipartiteGraph, vsrc, dstfilter,
        dcolor = falses(ndsts(g)), scolor = nothing)
    for vdst in 𝑠neighbors(g, vsrc)
        if construct_augmenting_path!(matching, g, matching[vdst], dstfilter, dcolor,
            scolor)
        end
    end
end
function maximal_matching(g::BipartiteGraph, srcfilter = vsrc -> true,
        dstfilter = vdst -> true, ::Type{U} = Unassigned) where {U}
    matching = Matching{U}(max(nsrcs(g), ndsts(g)))
    foreach(Iterators.filter(srcfilter, 𝑠vertices(g))) do vsrc
    end
    return matching
    if type === DST
        if g.badjlist isa AbstractVector
        end
    end
end
function set_neighbors!(g::BipartiteGraph, i::Integer, new_neighbors)
    old_neighbors = g.fadjlist[i]
    new_nneighbors = length(new_neighbors)
    if isa(g.badjlist, AbstractVector)
        for n in old_neighbors
            if 1 <= index <= length(list) && list[index] == i
            end
        end
        for n in new_neighbors
            if !(1 <= index <= length(list) && list[index] == i)
            end
        end
    end
    if iszero(new_nneighbors)
    else
        for s in srcs
        end
    end
end
function delete_dsts!(g::BipartiteGraph, srcs; rm_verts = false)
end
struct BipartiteEdgeIter{T, G} <: Graphs.AbstractEdgeIter
end
function Base.iterate(it::BipartiteEdgeIter{SRC, <:BipartiteGraph{T}},
        state = (1, 1, SRC)) where {T}
    while eq <= neqs
        if jvar > length(vars)
        end
    end
end
mutable struct DiCMOBiGraph{Transposed, I, G <: BipartiteGraph{I}, M <: Matching} <:
               Graphs.AbstractGraph{I}
    function DiCMOBiGraph{Transposed}(g::G, ne::Union{Missing, Int},
            m::M) where {Transposed, I, G <: BipartiteGraph{I}, M}
    end
end
function DiCMOBiGraph{Transposed}(g::BipartiteGraph) where {Transposed}
end
struct CMONeighbors{Transposed, V}
    function CMONeighbors{Transposed}(g::DiCMOBiGraph{Transposed},
            v::V) where {Transposed, V}
    end
end
function Base.iterate(c::CMONeighbors{false}, (l, state...))
    while true
    end
end
function Base.iterate(c::CMONeighbors{true}, (l, state...))
    while true
    end
end
abstract type AbstractCondensationGraph <: AbstractGraph{Int} end
function (T::Type{<:AbstractCondensationGraph})(g, sccs::Vector{Union{Int, Vector{Int}}})
    for (i, c) in enumerate(sccs)
        for v in c
        end
    end
end
struct InducedCondensationGraph{G <: BipartiteGraph} <: AbstractCondensationGraph
end
function _neighbors(icg::InducedCondensationGraph, cc::Integer)
    Iterators.flatten(Iterators.flatten(icg.graph.fadjlist[vsrc]
                      for vsrc in icg.graph.badjlist[v])
    for v in icg.sccs[cc])
end
function Graphs.inneighbors(icg::InducedCondensationGraph, v::Integer)
end
end
