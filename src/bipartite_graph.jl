module BipartiteGraphs
import ModelingToolkit: complete
export BipartiteEdge, BipartiteGraph, DiCMOBiGraph, Unassigned, unassigned,
       Matching, InducedCondensationGraph, maximal_matching,
       construct_augmenting_path!, MatchedCondensationGraph
export 𝑠vertices, 𝑑vertices, has_𝑠vertex, has_𝑑vertex, 𝑠neighbors, 𝑑neighbors,
       𝑠edges, 𝑑edges, nsrcs, ndsts, SRC, DST, set_neighbors!, invview,
       delete_srcs!, delete_dsts!
using DocStringExtensions
using UnPack
using SparseArrays
using Graphs
using Setfield
struct Unassigned
    global unassigned
    const unassigned = Unassigned.instance
end
Base.length(u::Unassigned) = 1
Base.size(u::Unassigned) = ()
Base.iterate(u::Unassigned) = (unassigned, nothing)
Base.iterate(u::Unassigned, state) = nothing
Base.show(io::IO, ::Unassigned) = printstyled(io, "u"; color = :light_black)
struct Matching{U, V <: AbstractVector} <: AbstractVector{Union{U, Int}}
    match::V
    inv_match::Union{Nothing, V}
end
function Matching{V}(m::Matching) where {V}
    eltype(m) === Union{V, Int} && return M
    VUT = typeof(similar(m.match, Union{V, Int}))
    Matching{V}(convert(VUT, m.match),
        m.inv_match === nothing ? nothing : convert(VUT, m.inv_match))
end
Matching(m::Matching) = m
Matching{U}(v::V) where {U, V <: AbstractVector} = Matching{U, V}(v, nothing)
function Matching{U}(v::V, iv::Union{V, Nothing}) where {U, V <: AbstractVector}
    Matching{U, V}(v, iv)
end
function Matching(v::V) where {U, V <: AbstractVector{Union{U, Int}}}
    Matching{@isdefined(U) ? U : Unassigned, V}(v, nothing)
end
function Matching(m::Int)
    Matching{Unassigned}(Union{Int, Unassigned}[unassigned for _ in 1:m], nothing)
end
function Matching{U}(m::Int) where {U}
    Matching{Union{Unassigned, U}}(Union{Int, Unassigned, U}[unassigned for _ in 1:m],
        nothing)
end
Base.size(m::Matching) = Base.size(m.match)
Base.getindex(m::Matching, i::Integer) = m.match[i]
Base.iterate(m::Matching, state...) = iterate(m.match, state...)
function Base.copy(m::Matching{U}) where {U}
    Matching{U}(copy(m.match), m.inv_match === nothing ? nothing : copy(m.inv_match))
end
function Base.setindex!(m::Matching{U}, v::Union{Integer, U}, i::Integer) where {U}
    if m.inv_match !== nothing
        oldv = m.match[i]
        if v isa Int && 1 <= v <= length(m.inv_match) && (iv = m.inv_match[v]) isa Int
            m.match[iv] = unassigned
        end
        if isa(oldv, Int)
            @assert m.inv_match[oldv] == i
            m.inv_match[oldv] = unassigned
        end
        if isa(v, Int)
            for vv in (length(m.inv_match) + 1):v
                push!(m.inv_match, unassigned)
            end
            m.inv_match[v] = i
        end
    end
    return m.match[i] = v
end
function Base.push!(m::Matching, v)
    push!(m.match, v)
    if v isa Integer && m.inv_match !== nothing
        for vv in (length(m.inv_match) + 1):v
            push!(m.inv_match, unassigned)
        end
        m.inv_match[v] = length(m.match)
    end
end
function complete(m::Matching{U},
        N = maximum((x for x in m.match if isa(x, Int)); init = 0)) where {U}
    m.inv_match !== nothing && return m
    inv_match = Union{U, Int}[unassigned for _ in 1:N]
    for (i, eq) in enumerate(m.match)
        isa(eq, Int) || continue
        inv_match[eq] = i
    end
    return Matching{U}(collect(m.match), inv_match)
end
@noinline function require_complete(m::Matching)
    m.inv_match === nothing &&
        throw(ArgumentError("Backwards matching not defined. `complete` the matching first."))
end
function invview(m::Matching{U, V}) where {U, V}
    require_complete(m)
    return Matching{U, V}(m.inv_match, m.match)
end
@enum VertType SRC DST
struct BipartiteEdge{I <: Integer} <: Graphs.AbstractEdge{I}
    src::I
    dst::I
    function BipartiteEdge(src::I, dst::V) where {I, V}
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
""""""
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
        badj::Union{AbstractVector, Integer} = maximum(maximum, fadj);
        metadata = nothing)
    BipartiteGraph(mapreduce(length, +, fadj; init = 0), fadj, badj, metadata)
end
@noinline function require_complete(g::BipartiteGraph)
    g.badjlist isa AbstractVector ||
        throw(ArgumentError("The graph has no back edges. Use `complete`."))
end
function invview(g::BipartiteGraph)
    require_complete(g)
    BipartiteGraph(g.ne, g.badjlist, g.fadjlist)
end
function complete(g::BipartiteGraph{I}) where {I}
    isa(g.badjlist, AbstractVector) && return g
    badjlist = Vector{I}[Vector{I}() for _ in 1:(g.badjlist)]
    for (s, l) in enumerate(g.fadjlist)
        for d in l
            push!(badjlist[d], s)
        end
    end
    BipartiteGraph(g.ne, g.fadjlist, badjlist)
end
struct BipartiteAdjacencyList
    u::Union{Vector{Int}, Nothing}
    highlight_u::Union{Set{Int}, Nothing}
    match::Union{Int, Bool, Unassigned}
end
function BipartiteAdjacencyList(u::Union{Vector{Int}, Nothing})
    BipartiteAdjacencyList(u, nothing, unassigned)
end
struct HighlightInt
    i::Int
    highlight::Symbol
    match::Bool
end
Base.typeinfo_implicit(::Type{HighlightInt}) = true
function Base.show(io::IO, hi::HighlightInt)
    if hi.match
        printstyled(io, "(", color = hi.highlight)
        printstyled(io, hi.i, color = hi.highlight)
        printstyled(io, ")", color = hi.highlight)
    else
        printstyled(io, hi.i, color = hi.highlight)
    end
end
function Base.show(io::IO, l::BipartiteAdjacencyList)
    if l.match === true
        printstyled(io, "∫ ", color = :cyan)
    else
        printstyled(io, "  ")
    end
    if l.u === nothing
        printstyled(io, '⋅', color = :light_black)
    elseif isempty(l.u)
        printstyled(io, '∅', color = :light_black)
    elseif l.highlight_u === nothing
        print(io, l.u)
    else
        match = l.match
        isa(match, Bool) && (match = unassigned)
        function choose_color(i)
            solvable = i in l.highlight_u
            matched = i == match
            if !matched && solvable
                :default
            elseif !matched && !solvable
                :light_black
            elseif matched && solvable
                :light_yellow
            elseif matched && !solvable
                :magenta
            end
        end
        if !isempty(setdiff(l.highlight_u, l.u))
            print(io,
                map(union(l.u, l.highlight_u)) do i
                    HighlightInt(i, !(i in l.u) ? :light_red : choose_color(i),
                        i == match)
                end)
        else
            print(io, map(l.u) do i
                HighlightInt(i, choose_color(i), i == match)
            end)
        end
    end
end
struct Label
    s::String
    c::Symbol
end
Label(s::AbstractString) = Label(s, :nothing)
Label(x::Integer) = Label(string(x))
Base.show(io::IO, l::Label) = printstyled(io, l.s, color = l.c)
struct BipartiteGraphPrintMatrix <:
       AbstractMatrix{Union{Label, Int, BipartiteAdjacencyList}}
    bpg::BipartiteGraph
end
Base.size(bgpm::BipartiteGraphPrintMatrix) = (max(nsrcs(bgpm.bpg), ndsts(bgpm.bpg)) + 1, 3)
function Base.getindex(bgpm::BipartiteGraphPrintMatrix, i::Integer, j::Integer)
    checkbounds(bgpm, i, j)
    if i == 1
        return (Label.(("#", "src", "dst")))[j]
    elseif j == 1
        return i - 1
    elseif j == 2
        return BipartiteAdjacencyList(i - 1 <= nsrcs(bgpm.bpg) ?
                                      𝑠neighbors(bgpm.bpg, i - 1) : nothing)
    elseif j == 3
        return BipartiteAdjacencyList(i - 1 <= ndsts(bgpm.bpg) ?
                                      𝑑neighbors(bgpm.bpg, i - 1) : nothing)
    else
        @assert false
    end
end
function Base.show(io::IO, b::BipartiteGraph)
    print(io, "BipartiteGraph with (", length(b.fadjlist), ", ",
        isa(b.badjlist, Int) ? b.badjlist : length(b.badjlist), ") (𝑠,𝑑)-vertices\n")
    Base.print_matrix(io, BipartiteGraphPrintMatrix(b))
end
""""""
function Base.isequal(bg1::BipartiteGraph{T}, bg2::BipartiteGraph{T}) where {T <: Integer}
    iseq = (bg1.ne == bg2.ne)
    iseq &= (bg1.fadjlist == bg2.fadjlist)
    iseq &= (bg1.badjlist == bg2.badjlist)
    iseq
end
""""""
function BipartiteGraph(nsrcs::T, ndsts::T, backedge::Val{B} = Val(true);
        metadata = nothing) where {T, B}
    fadjlist = map(_ -> T[], 1:nsrcs)
    badjlist = B ? map(_ -> T[], 1:ndsts) : ndsts
    BipartiteGraph(0, fadjlist, badjlist, metadata)
end
function Base.copy(bg::BipartiteGraph)
    BipartiteGraph(bg.ne, map(copy, bg.fadjlist), map(copy, bg.badjlist),
        deepcopy(bg.metadata))
end
Base.eltype(::Type{<:BipartiteGraph{I}}) where {I} = I
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
if isdefined(Graphs, :has_contiguous_vertices)
    Graphs.has_contiguous_vertices(::Type{<:BipartiteGraph}) = false
end
Graphs.is_directed(::Type{<:BipartiteGraph}) = false
Graphs.vertices(g::BipartiteGraph) = (𝑠vertices(g), 𝑑vertices(g))
𝑠vertices(g::BipartiteGraph) = axes(g.fadjlist, 1)
function 𝑑vertices(g::BipartiteGraph)
    g.badjlist isa AbstractVector ? axes(g.badjlist, 1) : Base.OneTo(g.badjlist)
end
has_𝑠vertex(g::BipartiteGraph, v::Integer) = v in 𝑠vertices(g)
has_𝑑vertex(g::BipartiteGraph, v::Integer) = v in 𝑑vertices(g)
function 𝑠neighbors(g::BipartiteGraph, i::Integer,
        with_metadata::Val{M} = Val(false)) where {M}
    M ? zip(g.fadjlist[i], g.metadata[i]) : g.fadjlist[i]
end
function 𝑑neighbors(g::BipartiteGraph, j::Integer,
        with_metadata::Val{M} = Val(false)) where {M}
    require_complete(g)
    backj = g.badjlist[j]::Vector{Int}
    M ? zip(backj, (g.metadata[i][j] for i in backj)) : backj
end
Graphs.ne(g::BipartiteGraph) = g.ne
Graphs.nv(g::BipartiteGraph) = sum(length, vertices(g))
Graphs.edgetype(g::BipartiteGraph{I}) where {I} = BipartiteEdge{I}
nsrcs(g::BipartiteGraph) = length(𝑠vertices(g))
ndsts(g::BipartiteGraph) = length(𝑑vertices(g))
function Graphs.has_edge(g::BipartiteGraph, edge::BipartiteEdge)
    @unpack src, dst = edge
    (src in 𝑠vertices(g) && dst in 𝑑vertices(g)) || return false
    insorted(dst, 𝑠neighbors(g, src))
end
Base.in(edge::BipartiteEdge, g::BipartiteGraph) = Graphs.has_edge(g, edge)
""""""
function construct_augmenting_path!(matching::Matching, g::BipartiteGraph, vsrc, dstfilter,
        dcolor = falses(ndsts(g)), scolor = nothing)
    scolor === nothing || (scolor[vsrc] = true)
    for vdst in 𝑠neighbors(g, vsrc)
        if dstfilter(vdst) && matching[vdst] === unassigned
            matching[vdst] = vsrc
            return true
        end
    end
    for vdst in 𝑠neighbors(g, vsrc)
        (dstfilter(vdst) && !dcolor[vdst]) || continue
        dcolor[vdst] = true
        if construct_augmenting_path!(matching, g, matching[vdst], dstfilter, dcolor,
            scolor)
            matching[vdst] = vsrc
            return true
        end
    end
    return false
end
""""""
function maximal_matching(g::BipartiteGraph, srcfilter = vsrc -> true,
        dstfilter = vdst -> true, ::Type{U} = Unassigned) where {U}
    matching = Matching{U}(max(nsrcs(g), ndsts(g)))
    foreach(Iterators.filter(srcfilter, 𝑠vertices(g))) do vsrc
        construct_augmenting_path!(matching, g, vsrc, dstfilter)
    end
    return matching
end
struct NoMetadata end
const NO_METADATA = NoMetadata()
function Graphs.add_edge!(g::BipartiteGraph, i::Integer, j::Integer, md = NO_METADATA)
    add_edge!(g, BipartiteEdge(i, j), md)
end
function Graphs.add_edge!(g::BipartiteGraph, edge::BipartiteEdge, md = NO_METADATA)
    @unpack fadjlist, badjlist = g
    s, d = src(edge), dst(edge)
    (has_𝑠vertex(g, s) && has_𝑑vertex(g, d)) || error("edge ($edge) out of range.")
    @inbounds list = fadjlist[s]
    index = searchsortedfirst(list, d)
    @inbounds (index <= length(list) && list[index] == d) && return false
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
    return true
end
function Graphs.rem_edge!(g::BipartiteGraph, i::Integer, j::Integer)
    Graphs.rem_edge!(g, BipartiteEdge(i, j))
end
function Graphs.rem_edge!(g::BipartiteGraph, edge::BipartiteEdge)
    @unpack fadjlist, badjlist = g
    s, d = src(edge), dst(edge)
    (has_𝑠vertex(g, s) && has_𝑑vertex(g, d)) || error("edge ($edge) out of range.")
    @inbounds list = fadjlist[s]
    index = searchsortedfirst(list, d)
    @inbounds (index <= length(list) && list[index] == d) ||
              error("graph does not have edge $edge")
    deleteat!(list, index)
    g.ne -= 1
    if badjlist isa AbstractVector
        @inbounds list = badjlist[d]
        index = searchsortedfirst(list, s)
        deleteat!(list, index)
    end
    return true
end
function Graphs.add_vertex!(g::BipartiteGraph{T}, type::VertType) where {T}
    if type === DST
        if g.badjlist isa AbstractVector
            push!(g.badjlist, T[])
            return length(g.badjlist)
        else
            g.badjlist += 1
            return g.badjlist
        end
    elseif type === SRC
        push!(g.fadjlist, T[])
        return length(g.fadjlist)
    else
        error("type ($type) must be either `DST` or `SRC`")
    end
end
function set_neighbors!(g::BipartiteGraph, i::Integer, new_neighbors)
    old_neighbors = g.fadjlist[i]
    old_nneighbors = length(old_neighbors)
    new_nneighbors = length(new_neighbors)
    g.ne += new_nneighbors - old_nneighbors
    if isa(g.badjlist, AbstractVector)
        for n in old_neighbors
            @inbounds list = g.badjlist[n]
            index = searchsortedfirst(list, i)
            if 1 <= index <= length(list) && list[index] == i
                deleteat!(list, index)
            end
        end
        for n in new_neighbors
            @inbounds list = g.badjlist[n]
            index = searchsortedfirst(list, i)
            if !(1 <= index <= length(list) && list[index] == i)
                insert!(list, index, i)
            end
        end
    end
    if iszero(new_nneighbors)
        empty!(g.fadjlist[i])
    else
        g.fadjlist[i] = unique!(sort(new_neighbors))
    end
end
function delete_srcs!(g::BipartiteGraph{I}, srcs; rm_verts = false) where {I}
    for s in srcs
        set_neighbors!(g, s, ())
    end
    if rm_verts
        old_to_new_idxs = collect(one(I):I(nsrcs(g)))
        for s in srcs
            old_to_new_idxs[s] = zero(I)
        end
        offset = zero(I)
        for i in eachindex(old_to_new_idxs)
            if iszero(old_to_new_idxs[i])
                offset += one(I)
                continue
            end
            old_to_new_idxs[i] -= offset
        end
        if g.badjlist isa AbstractVector
            for i in 1:ndsts(g)
                for j in eachindex(g.badjlist[i])
                    g.badjlist[i][j] = old_to_new_idxs[g.badjlist[i][j]]
                end
                filter!(!iszero, g.badjlist[i])
            end
        end
        deleteat!(g.fadjlist, srcs)
    end
    g
end
function delete_dsts!(g::BipartiteGraph, srcs; rm_verts = false)
    delete_srcs!(invview(g), srcs; rm_verts)
end
Graphs.edges(g::BipartiteGraph) = BipartiteEdgeIter(g, Val(SRC))
𝑠edges(g::BipartiteGraph) = BipartiteEdgeIter(g, Val(SRC))
𝑑edges(g::BipartiteGraph) = BipartiteEdgeIter(g, Val(DST))
struct BipartiteEdgeIter{T, G} <: Graphs.AbstractEdgeIter
    g::G
    type::Val{T}
end
Base.length(it::BipartiteEdgeIter) = ne(it.g)
Base.eltype(it::BipartiteEdgeIter) = edgetype(it.g)
function Base.iterate(it::BipartiteEdgeIter{SRC, <:BipartiteGraph{T}},
        state = (1, 1, SRC)) where {T}
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
function Base.iterate(it::BipartiteEdgeIter{DST, <:BipartiteGraph{T}},
        state = (1, 1, DST)) where {T}
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
function Graphs.incidence_matrix(g::BipartiteGraph, val = true)
    I = Int[]
    J = Int[]
    for i in 𝑠vertices(g), n in 𝑠neighbors(g, i)
        push!(I, i)
        push!(J, n)
    end
    S = sparse(I, J, val, nsrcs(g), ndsts(g))
end
""""""
mutable struct DiCMOBiGraph{Transposed, I, G <: BipartiteGraph{I}, M <: Matching} <:
               Graphs.AbstractGraph{I}
    graph::G
    ne::Union{Missing, Int}
    matching::M
    function DiCMOBiGraph{Transposed}(g::G, ne::Union{Missing, Int},
            m::M) where {Transposed, I, G <: BipartiteGraph{I}, M}
        new{Transposed, I, G, M}(g, ne, m)
    end
end
function DiCMOBiGraph{Transposed}(g::BipartiteGraph) where {Transposed}
    DiCMOBiGraph{Transposed}(g, 0, Matching(ndsts(g)))
end
function DiCMOBiGraph{Transposed}(g::BipartiteGraph, m::M) where {Transposed, M}
    DiCMOBiGraph{Transposed}(g, missing, m)
end
function invview(g::DiCMOBiGraph{Transposed}) where {Transposed}
    DiCMOBiGraph{!Transposed}(invview(g.graph), g.ne, invview(g.matching))
end
Graphs.is_directed(::Type{<:DiCMOBiGraph}) = true
function Graphs.nv(g::DiCMOBiGraph{Transposed}) where {Transposed}
    Transposed ? ndsts(g.graph) : nsrcs(g.graph)
end
function Graphs.vertices(g::DiCMOBiGraph{Transposed}) where {Transposed}
    Transposed ? 𝑑vertices(g.graph) : 𝑠vertices(g.graph)
end
struct CMONeighbors{Transposed, V}
    g::DiCMOBiGraph{Transposed}
    v::V
    function CMONeighbors{Transposed}(g::DiCMOBiGraph{Transposed},
            v::V) where {Transposed, V}
        new{Transposed, V}(g, v)
    end
end
Graphs.outneighbors(g::DiCMOBiGraph{false}, v) = CMONeighbors{false}(g, v)
Graphs.inneighbors(g::DiCMOBiGraph{false}, v) = inneighbors(invview(g), v)
Base.iterate(c::CMONeighbors{false}) = iterate(c, (c.g.graph.fadjlist[c.v],))
function Base.iterate(c::CMONeighbors{false}, (l, state...))
    while true
        r = iterate(l, state...)
        r === nothing && return nothing
        vsrc = c.g.matching[r[1]]
        if vsrc === c.v || !isa(vsrc, Int)
            state = (r[2],)
            continue
        end
        return vsrc, (l, r[2])
    end
end
Base.length(c::CMONeighbors{false}) = count(_ -> true, c)
liftint(f, x) = (!isa(x, Int)) ? nothing : f(x)
liftnothing(f, x) = x === nothing ? nothing : f(x)
_vsrc(c::CMONeighbors{true}) = c.g.matching[c.v]
_neighbors(c::CMONeighbors{true}) = liftint(vsrc -> c.g.graph.fadjlist[vsrc], _vsrc(c))
Base.length(c::CMONeighbors{true}) = something(liftnothing(length, _neighbors(c)), 1) - 1
Graphs.inneighbors(g::DiCMOBiGraph{true}, v) = CMONeighbors{true}(g, v)
Graphs.outneighbors(g::DiCMOBiGraph{true}, v) = outneighbors(invview(g), v)
Base.iterate(c::CMONeighbors{true}) = liftnothing(ns -> iterate(c, (ns,)), _neighbors(c))
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
function _edges(g::DiCMOBiGraph{Transposed}) where {Transposed}
    Transposed ?
    ((w => v for w in inneighbors(g, v)) for v in vertices(g)) :
    ((v => w for w in outneighbors(g, v)) for v in vertices(g))
end
Graphs.edges(g::DiCMOBiGraph) = (Graphs.SimpleEdge(p) for p in Iterators.flatten(_edges(g)))
function Graphs.ne(g::DiCMOBiGraph)
    if g.ne === missing
        g.ne = mapreduce(x -> length(x.iter), +, _edges(g))
    end
    return g.ne
end
Graphs.has_edge(g::DiCMOBiGraph{true}, a, b) = a in inneighbors(g, b)
Graphs.has_edge(g::DiCMOBiGraph{false}, a, b) = b in outneighbors(g, a)
(::Type{<:DiCMOBiGraph})(n::Integer) = SimpleDiGraph(n)
abstract type AbstractCondensationGraph <: AbstractGraph{Int} end
function (T::Type{<:AbstractCondensationGraph})(g, sccs::Vector{Union{Int, Vector{Int}}})
    scc_assignment = Vector{Int}(undef, isa(g, BipartiteGraph) ? ndsts(g) : nv(g))
    for (i, c) in enumerate(sccs)
        for v in c
            scc_assignment[v] = i
        end
    end
    T(g, sccs, scc_assignment)
end
function (T::Type{<:AbstractCondensationGraph})(g, sccs::Vector{Vector{Int}})
    T(g, Vector{Union{Int, Vector{Int}}}(sccs))
end
Graphs.is_directed(::Type{<:AbstractCondensationGraph}) = true
Graphs.nv(icg::AbstractCondensationGraph) = length(icg.sccs)
Graphs.vertices(icg::AbstractCondensationGraph) = Base.OneTo(nv(icg))
""""""
struct MatchedCondensationGraph{G <: DiCMOBiGraph} <: AbstractCondensationGraph
    graph::G
    sccs::Vector{Union{Int, Vector{Int}}}
    scc_assignment::Vector{Int}
end
function Graphs.outneighbors(mcg::MatchedCondensationGraph, cc::Integer)
    Iterators.flatten((mcg.scc_assignment[v′]
                      for v′ in outneighbors(mcg.graph, v) if mcg.scc_assignment[v′] != cc)
    for v in mcg.sccs[cc])
end
function Graphs.inneighbors(mcg::MatchedCondensationGraph, cc::Integer)
    Iterators.flatten((mcg.scc_assignment[v′]
                      for v′ in inneighbors(mcg.graph, v) if mcg.scc_assignment[v′] != cc)
    for v in mcg.sccs[cc])
end
""""""
struct InducedCondensationGraph{G <: BipartiteGraph} <: AbstractCondensationGraph
    graph::G
    sccs::Vector{Union{Int, Vector{Int}}}
    scc_assignment::Vector{Int}
end
function _neighbors(icg::InducedCondensationGraph, cc::Integer)
    Iterators.flatten(Iterators.flatten(icg.graph.fadjlist[vsrc]
                      for vsrc in icg.graph.badjlist[v])
    for v in icg.sccs[cc])
end
function Graphs.outneighbors(icg::InducedCondensationGraph, v::Integer)
    (icg.scc_assignment[n] for n in _neighbors(icg, v) if icg.scc_assignment[n] > v)
end
function Graphs.inneighbors(icg::InducedCondensationGraph, v::Integer)
    (icg.scc_assignment[n] for n in _neighbors(icg, v) if icg.scc_assignment[n] < v)
end
end
