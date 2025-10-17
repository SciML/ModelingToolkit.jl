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
function Matching{U}(m::Int) where {U}
    Matching{Union{Unassigned, U}}(Union{Int, Unassigned, U}[unassigned for _ in 1:m],
        nothing)
end
Base.size(m::Matching) = Base.size(m.match)
Base.getindex(m::Matching, i::Integer) = m.match[i]
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
function complete(g::BipartiteGraph{I}) where {I}
    badjlist = Vector{I}[Vector{I}() for _ in 1:(g.badjlist)]
    BipartiteGraph(g.ne, g.fadjlist, badjlist)
end
function BipartiteGraph(nsrcs::T, ndsts::T, backedge::Val{B} = Val(true);
        metadata = nothing) where {T, B}
    fadjlist = map(_ -> T[], 1:nsrcs)
    badjlist = B ? map(_ -> T[], 1:ndsts) : ndsts
    BipartiteGraph(0, fadjlist, badjlist, metadata)
end
𝑠vertices(g::BipartiteGraph) = axes(g.fadjlist, 1)
function 𝑑vertices(g::BipartiteGraph)
    g.badjlist isa AbstractVector ? axes(g.badjlist, 1) : Base.OneTo(g.badjlist)
end
function 𝑠neighbors(g::BipartiteGraph, i::Integer,
        with_metadata::Val{M} = Val(false)) where {M}
    M ? zip(g.fadjlist[i], g.metadata[i]) : g.fadjlist[i]
end
nsrcs(g::BipartiteGraph) = length(𝑠vertices(g))
ndsts(g::BipartiteGraph) = length(𝑑vertices(g))
function maximal_matching(g::BipartiteGraph, srcfilter = vsrc -> true,
        dstfilter = vdst -> true, ::Type{U} = Unassigned) where {U}
    matching = Matching{U}(max(nsrcs(g), ndsts(g)))
    return matching
end
mutable struct DiCMOBiGraph{Transposed, I, G <: BipartiteGraph{I}, M <: Matching} <:
               Graphs.AbstractGraph{I}
end
end
