"""
    $(TYPEDEF)

A vertex in the connection hypergraph.

## Fields

$(TYPEDFIELDS)
"""
struct ConnectionVertex
    """
    The name of the variable or subsystem represented by this connection vertex. Stored as
    a list of names denoting the path from the root system to this variable/subsystem. The
    name of the root system is not included.
    """
    name::Vector{Symbol}
    """
    Boolean indicating whether this is an outside connector.
    """
    isouter::Bool
    """
    A type indicating what kind of connector it is. One of:
    - `Stream`
    - `Equality`
    - `Flow`
    - `InputVar`
    - `OutputVar`
    """
    type::DataType
    """
    The cached hash value of this struct. Should never be passed manually.
    """
    hash::UInt
end

"""
    $(TYPEDSIGNATURES)

Create a `ConnectionVertex` given
- `namespace`: the path from the root to the variable/subsystem. Does not include the root
  system.
- `var`: the variable/subsystem.

`isouter` is the same as the struct field. Uses `get_connection_type` to find the type to
use for this connection.
"""
function ConnectionVertex(
        namespace::Vector{Symbol}, var::Union{BasicSymbolic, AbstractSystem}, isouter::Bool)
    if var isa BasicSymbolic
        name = getname(var)
    else
        name = nameof(var)
    end
    var_ns = namespace_hierarchy(name)
    type = get_connection_type(var)
    name = vcat(namespace, var_ns)
    return ConnectionVertex(name, isouter, type; alias = true)
end

"""
    $(TYPEDSIGNATURES)

Create a connection vertex for the given path. Typically used for domain connection graphs,
where the type of connection doesn't matter. Uses `isouter = true` and `type = Flow`.
"""
function ConnectionVertex(name::Vector{Symbol})
    return ConnectionVertex(name, true, Flow)
end

"""
    $(TYPEDSIGNATURES)

Create a connection vertex for the given path `name` using the provided `isouter` and
`type`. `alias` denotes whether `name` can be stored by this vertex without copying.
"""
function ConnectionVertex(
        name::Vector{Symbol}, isouter::Bool, type::DataType; alias = false)
    if !alias
        name = copy(name)
    end
    h = foldr(hash, name; init = zero(UInt))
    h = hash(type, h)
    h = hash(isouter, h)
    return ConnectionVertex(name, isouter, type, h)
end

Base.hash(x::ConnectionVertex, h::UInt) = h ⊻ x.hash

function Base.:(==)(a::ConnectionVertex, b::ConnectionVertex)
    length(a.name) == length(b.name) || return false
    for (x, y) in zip(a.name, b.name)
        x == y || return false
    end
    a.isouter == b.isouter || return false
    a.type == b.type || return false
    if a.hash != b.hash
        error("""
        This should never happen. Please open an issue in ModelingToolkit.jl with an MWE.
        """)
    end
    return true
end

function Base.show(io::IO, vert::ConnectionVertex)
    for name in @view(vert.name[1:(end - 1)])
        print(io, name, ".")
    end
    print(io, vert.name[end], "::", vert.isouter ? "outer" : "inner")
end

"""
    $(TYPEDEF)

A hypergraph used to represent the connection sets in a system. Vertices of this graph are
of type `ConnectionVertex`. The connected components of a connection graph are the merged
connection sets.

## Fields

$(TYPEDFIELDS)
"""
struct HyperGraph{V}
    """
    Mapping from vertices to their integer ID.
    """
    labels::Dict{V, Int}
    """
    Reverse mapping from integer ID to vertices.
    """
    invmap::Vector{V}
    """
    Core data structure for storing the hypergraph. Each hyperedge is a source vertex and
    has bipartite edges to the connection vertices it is incident on.
    """
    graph::BipartiteGraph{Int, Nothing}
end

const ConnectionGraph = HyperGraph{ConnectionVertex}

"""
    $(TYPEDSIGNATURES)

Create an empty `ConnectionGraph`.
"""
function HyperGraph{V}() where {V}
    graph = BipartiteGraph(0, 0, Val(true))
    return HyperGraph{V}(Dict{V, Int}(), V[], graph)
end

function Base.show(io::IO, graph::ConnectionGraph)
    printstyled(io, get(io, :cgraph_name, "ConnectionGraph"); color = :blue, bold = true)
    println(io, " with ", length(graph.labels),
        " vertices and ", nsrcs(graph.graph), " hyperedges")
    compact = get(io, :compact, false)
    for edge_i in 𝑠vertices(graph.graph)
        if compact && edge_i > 5
            println(io, "⋮")
            break
        end
        edge_idxs = 𝑠neighbors(graph.graph, edge_i)
        type = graph.invmap[edge_idxs[1]].type
        if type <: Union{InputVar, OutputVar}
            type = "Causal"
        elseif type == Equality
            # otherwise it prints `ModelingToolkit.Equality`
            type = "Equality"
        end
        printstyled(io, "  ", type; bold = true, color = :yellow)
        print(io, "<")
        for vi in @view(edge_idxs[1:(end - 1)])
            print(io, graph.invmap[vi], ", ")
        end
        println(io, graph.invmap[edge_idxs[end]], ">")
    end
end

"""
    $(TYPEDSIGNATURES)

Add the given vertex to the connection graph. Return the integer ID of the added vertex.
No-op if the vertex already exists.
"""
function Graphs.add_vertex!(graph::HyperGraph{V}, dst::V) where {V}
    j = get(graph.labels, dst, 0)
    iszero(j) || return j
    j = Graphs.add_vertex!(graph.graph, DST)
    push!(graph.invmap, dst)
    @assert length(graph.invmap) == j
    graph.labels[dst] = j
    return j
end

const HyperGraphEdge{V} = Union{Vector{V}, Tuple{Vararg{V}}, Set{V}}
const ConnectionGraphEdge = HyperGraphEdge{ConnectionVertex}

"""
    $(TYPEDSIGNATURES)

Add the given hyperedge to the connection graph. Adds all vertices in the given edge if
they do not exist. Returns the integer ID of the added edge.
"""
function Graphs.add_edge!(graph::HyperGraph{V}, src::HyperGraphEdge{V}) where {V}
    i = Graphs.add_vertex!(graph.graph, SRC)
    for vert in src
        j = Graphs.add_vertex!(graph, vert)
        Graphs.add_edge!(graph.graph, i, j)
    end
    return i
end

"""
    $(TYPEDEF)

A connection state is a combination of two `ConnectionGraph`s, one for the connection sets
and the other for the domain network. The domain network is a graph of connected
subsystems. The connected components of the domain network denote connected domains that
share properties.
"""
abstract type AbstractConnectionState end

"""
    $(TYPEDEF)

The most trivial `AbstractConnectionState`.

## Fields

$(TYPEDFIELDS)
"""
struct ConnectionState <: AbstractConnectionState
    """
    The connection graph for connection sets.
    """
    connection_graph::ConnectionGraph
    """
    The connection graph for the domain network.
    """
    domain_connection_graph::ConnectionGraph
end

"""
    $(TYPEDSIGNATURES)

Create an empty `ConnectionState` with empty graphs.
"""
ConnectionState() = ConnectionState(ConnectionGraph(), ConnectionGraph())

function Base.show(io::IO, state::AbstractConnectionState)
    printstyled(io, typeof(state); bold = true, color = :green)
    println(io, " comprising of")
    ctx1 = IOContext(io, :cgraph_name => "Connection Network", :compact => true)
    show(ctx1, state.connection_graph)
    println(io)
    println(io, "And")
    println(io)
    ctx2 = IOContext(io, :cgraph_name => "Domain Network", :compact => true)
    show(ctx2, state.domain_connection_graph)
end

"""
    $(TYPEDSIGNATURES)

Add the given edge to the connection network. Does not affect the domain network.
"""
function add_connection_edge!(state::ConnectionState, edge::ConnectionGraphEdge)
    Graphs.add_edge!(state.connection_graph, edge)
    return nothing
end

"""
    $(TYPEDSIGNATURES)

Add the given edge to the domain network. Does not affect the connection network.
"""
function add_domain_connection_edge!(state::ConnectionState, edge::ConnectionGraphEdge)
    Graphs.add_edge!(state.domain_connection_graph, edge)
    return nothing
end

"""
    $(TYPEDSIGNATURES)

An `AbstractConnectionState` that is used to remove edges from the main connection state.
Transformed analysis points add to the list of removed connections, and the list of removed
connections builds this connection state. This allows ensuring that the removed connections
are not present in the final network even if they are connected multiple times. This state
also tracks which vertex in each hyperedge is the input, since the removed connections are
causal.

## Fields

$(TYPEDFIELDS)
"""
struct NegativeConnectionState <: AbstractConnectionState
    """
    The connection graph for connection sets.
    """
    connection_graph::ConnectionGraph
    """
    The connection graph for the domain network.
    """
    domain_connection_graph::ConnectionGraph
    """
    Mapping from the integer ID of each hyperedge in `connection_graph` to the integer ID
    of the "input" in that hyperedge.
    """
    connection_hyperedge_inputs::Vector{Int}
    """
    Mapping from the integer ID of each hyperedge in `domain_connection_graph` to the
    integer ID of the "input" in that hyperedge.
    """
    domain_hyperedge_inputs::Vector{Int}
end

"""
    $(TYPEDSIGNATURES)

Create an empty `NegativeConnectionState` with empty graphs.
"""
function NegativeConnectionState()
    NegativeConnectionState(ConnectionGraph(), ConnectionGraph(), Int[], Int[])
end

"""
    $(TYPEDSIGNATURES)

Add the given edge to the connection network. Does not affect the domain network. Assumes
that the first vertex in `edge` is the input.
"""
function add_connection_edge!(state::NegativeConnectionState, edge::ConnectionGraphEdge)
    i = Graphs.add_edge!(state.connection_graph, edge)
    j = state.connection_graph.labels[first(edge)]
    push!(state.connection_hyperedge_inputs, j)
    @assert length(state.connection_hyperedge_inputs) == i
    return nothing
end

"""
    $(TYPEDSIGNATURES)

Add the given edge to the domain network. Does not affect the connection network. Assumes
that the first vertex in `edge` is the input.
"""
function add_domain_connection_edge!(
        state::NegativeConnectionState, edge::ConnectionGraphEdge)
    i = Graphs.add_edge!(state.domain_connection_graph, edge)
    j = state.domain_connection_graph.labels[first(edge)]
    push!(state.domain_hyperedge_inputs, j)
    @assert length(state.domain_hyperedge_inputs) == i
    return nothing
end

"""
    $(TYPEDSIGNATURES)

Modify `graph` such that no hyperedge is a superset of any (causal) hyerpedge in `neg_graph`.

For each "negative" hyperedge in `neg_graph` with integer ID `neg_edge_id`,
`neg_edge_inputs[neg_edge_id]` denotes the vertex the negative hyperedge is incident on
which is considered the input of the negative hyperedge. If any hyperedge in `graph`
contains this input as well as at least one other vertex in the negative hyperedge, all
vertices common between the hyperedge and negative hyperedge are removed from the hyperedge.

`graph` is modified in-place. Note that `graph` and `neg_graph` may not have the same
ordering of vertices, and thus all comparisons should be done by comparing the
`ConnectionVertex`.
"""
function remove_negative_connections!(
        graph::ConnectionGraph, neg_graph::ConnectionGraph, neg_edge_inputs::Vector{Int})
    # _i means index in neg_graph
    # _j means index in graph

    # get all edges in `graph` as bitsets
    graph_hyperedgesets = map(𝑠vertices(graph.graph)) do edge_j
        hyperedge_jdxs = 𝑠neighbors(graph.graph, edge_j)
        return BitSet(hyperedge_jdxs)
    end

    # indexes in each hyperedge to remove
    idxs_to_rm = [BitSet() for _ in graph_hyperedgesets]
    # iterate over negative edges and the corresponding input vertex in each edge
    for (input_i, edge_i) in zip(neg_edge_inputs, 𝑠vertices(neg_graph.graph))
        # get the hyperedge as integer indexes in `neg_graph`
        neg_hyperedge_idxs = 𝑠neighbors(neg_graph.graph, edge_i)
        # the hyperedge as `ConnectionVar`s
        neg_hyperedge = map(Base.Fix1(getindex, neg_graph.invmap), neg_hyperedge_idxs)
        # The hyperedge as integer indexes in `graph`
        # *j*dxs. See what I did there?
        neg_hyperedge_jdxs = map(cvar -> get(graph.labels, cvar, 0), neg_hyperedge)
        # the edge to remove is between variables that aren't connected, so ignore it
        if any(iszero, neg_hyperedge_jdxs)
            continue
        end

        # The input vertex as a `ConnectionVar`
        input_v = neg_graph.invmap[input_i]
        # The input vertex as an index in `graph`
        input_j = graph.labels[input_v]
        # Iterate over hyperedges in `graph`
        for edge_j in 𝑠vertices(graph.graph)
            # The bitset of nodes this edge is incident on
            edgeset = graph_hyperedgesets[edge_j]
            # the input must be in this hyperedge
            input_j in edgeset || continue
            # now, if any other vertex apart from this input is also in the hyperedge
            # we remove all the indices in `neg_hyperedge_jdxs` also present in this
            # hyperedge

            # should_rm tracks if any other vertex apart from `input_j` is in the hyperedge
            should_rm = false
            # iterate over the negative hyperedge
            for var_j in neg_hyperedge_jdxs
                var_j == input_j && continue
                # check if the variable which is not `input_j` is in the hyperedge
                should_rm |= var_j in edgeset
                should_rm || continue
                # if there is any other variable, start removing
                push!(idxs_to_rm[edge_j], var_j)
            end
        end
    end

    # for each edge and list of vertices to remove from the edge
    for (edge_j, neg_vertices) in enumerate(idxs_to_rm)
        for vert_j in neg_vertices
            # remove those vertices
            Graphs.rem_edge!(graph.graph, edge_j, vert_j)
        end
    end
end

"""
    $(TYPEDSIGNATURES)

Remove negative hyperedges given by `neg_state` from the connection and domain networks of
`state`.
"""
function remove_negative_connections!(
        state::ConnectionState, neg_state::NegativeConnectionState)
    remove_negative_connections!(state.connection_graph, neg_state.connection_graph,
        neg_state.connection_hyperedge_inputs)
    remove_negative_connections!(
        state.domain_connection_graph, neg_state.domain_connection_graph,
        neg_state.domain_hyperedge_inputs)
end

"""
    $(TYPEDSIGNATURES)

Return the merged connection sets in `graph` as a `Vector{Vector{ConnectionVertex}}`. These
are equivalent to the connected components of `graph`.
"""
function connectionsets(graph::HyperGraph{V}) where {V}
    bigraph = graph.graph
    invmap = graph.invmap

    # union all of the hyperedges
    disjoint_sets = IntDisjointSet(length(invmap))
    for edge_i in 𝑠vertices(bigraph)
        hyperedge = 𝑠neighbors(bigraph, edge_i)
        isempty(hyperedge) && continue
        root, rest = Iterators.peel(hyperedge)
        for vert in rest
            union!(disjoint_sets, root, vert)
        end
    end

    # maps the root of a vertex in `disjoint_sets` to the index of the corresponding set
    # in `vertex_sets`
    root_to_set = Dict{Int, Int}()
    vertex_sets = Vector{V}[]
    for (vert_i, vert) in enumerate(invmap)
        root = find_root!(disjoint_sets, vert_i)
        set_i = get!(root_to_set, root) do
            push!(vertex_sets, V[])
            return length(vertex_sets)
        end
        push!(vertex_sets[set_i], vert)
    end

    return vertex_sets
end

"""
    $(TYPEDSIGNATURES)

Return the connection sets of the connection graph and domain network.
"""
function connectionsets(state::ConnectionState)
    return connectionsets(state.connection_graph),
    connectionsets(state.domain_connection_graph)
end
