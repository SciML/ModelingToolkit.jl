module BipartiteGraphs
import ModelingToolkit: complete
export BipartiteGraph, Unassigned, unassigned,
       delete_srcs!, delete_dsts!
using Graphs
struct Unassigned
end
struct Matching{U, V <: AbstractVector} <: AbstractVector{Union}
end
function Matching(m::Matching) where {V}
    Matching(convert(VUT, m.match),
        m.inv_match === nothing ? nothing : convert(VUT, m.inv_match))
end
function Matching(v::V, iv::Union) where {U, V <: AbstractVector}
end
mutable struct BipartiteGraph{I <: Integer, M} <: Graphs.AbstractGraph{I}
    ne::Int
    fadjlist::Vector
    badjlist::Union{Vector{Vector{I}}, I}
    metadata::M
end
function BipartiteGraph(ne::Integer, fadj::AbstractVector,
        badj::Union{AbstractVector, Integer} = maximum(maximum, fadj);
        metadata = nothing)
    BipartiteGraph(ne, fadj, badj, metadata)
end
function complete(g::BipartiteGraph{I}) where {I}
    badjlist = Vector{I}[Vector() for _ in 1:(g.badjlist)]
    BipartiteGraph(g.ne, g.fadjlist, badjlist)
end
function BipartiteGraph(nsrcs::T, ndsts::T, backedge::Val{B} = Val(true);
        metadata = nothing) where {T, B}
    fadjlist = map(_ -> T[], 1:nsrcs)
    badjlist = B ? map(_ -> T[], 1:ndsts) : ndsts
    BipartiteGraph(0, fadjlist, badjlist, metadata)
end
function maximal_matching(g::BipartiteGraph, srcfilter = vsrc -> true,
        dstfilter = vdst -> true, ::Type = Unassigned) where {U}
end
end
