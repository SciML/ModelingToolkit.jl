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
mutable struct BipartiteGraph{I <: Integer, M} <: Graphs.AbstractGraph{I}
    ne::Int
    fadjlist::Vector
    badjlist::Union{Vector{Vector{I}}, I}
    metadata::M
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
