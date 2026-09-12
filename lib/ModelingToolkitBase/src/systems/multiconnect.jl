"""
    $(TYPEDEF)

A single labeled edge of a [`ConnectionNetwork`](@ref). The edge connects port
`src_port` of `nodes[src]` to port `dst_port` of `nodes[dst]`, where `nodes` is the
vector of node subsystems in the network. A port may be:

- a connector subsystem of the node (e.g. an acausal `Pin`), or
- an input/output variable of the node, for causal connections.

A port of `nothing` denotes the node itself; the node must be a connector in that
case. `src` and `dst` are 1-based indices into the `nodes` vector. For acausal
(undirected) connections the ordering of `src` and `dst` is immaterial. For causal
variable connections the convention is that `src` drives `dst`; the direction is in
any case determined by the input/output metadata of the connected variables.

See also [`multiconnect`](@ref), [`ConnectionNetwork`](@ref).
"""
struct ConnectionEdge
    src::Int
    dst::Int
    src_port::Union{Symbol, Nothing}
    dst_port::Union{Symbol, Nothing}
end

"""
    $(TYPEDSIGNATURES)

Create a `ConnectionEdge` between `src` and `dst` where the same port name `port` is
used on both ends.
"""
function ConnectionEdge(src::Integer, dst::Integer, port::Union{Symbol, Nothing})
    return ConnectionEdge(src, dst, port, port)
end

function Base.show(io::IO, edge::ConnectionEdge)
    _port_str(port) = port === nothing ? "self" : string(port)
    return print(
        io, "ConnectionEdge(", edge.src, " -> ", edge.dst, ", ",
        _port_str(edge.src_port), " -> ", _port_str(edge.dst_port), ")"
    )
end

"""
    $(TYPEDEF)

Payload of the `Connection` equation created by [`multiconnect`](@ref) when called
with a vector of node subsystems. Stores the instantiated node systems and the
labeled edges connecting their ports. The network is expanded into ordinary
connection sets by [`expand_connections`](@ref), where each edge behaves like a
pairwise `connect` between the two indicated ports.

# Fields

$(TYPEDFIELDS)
"""
struct ConnectionNetwork
    """
    The node subsystems taking part in the network. Each element is one
    instantiation of one of the (possibly several) component forms.
    """
    nodes::Vector{System}
    """
    The labeled edges between node ports.
    """
    edges::Vector{ConnectionEdge}
end

function Base.show(io::IO, net::ConnectionNetwork)
    return print(
        io, "ConnectionNetwork(", length(net.nodes), " nodes, ",
        length(net.edges), " edges)"
    )
end

"""
    $(TYPEDSIGNATURES)

Resolve port `port` of `node` to a connector subsystem or symbolic variable. A
`port` of `nothing` resolves to the node itself.
"""
function _network_port(node::System, port::Union{Symbol, Nothing})
    port === nothing && return node
    result = try
        getproperty(node, port)
    catch
        throw(
            ArgumentError(
            lazy"Node `$(nameof(node))` has no connector or variable named `$port`."
        )
        )
    end
    return result isa AbstractSystem ? result : unwrap(result)
end

"""
    $(TYPEDSIGNATURES)

Return the name of the port to use for `node` when the `portmap` of `multiconnect`
is `nothing`. If the node is itself a connector, returns `nothing`. Otherwise the
node must have exactly one connector subsystem, whose name is returned.
"""
function _default_network_port(node::System)
    connectors = Symbol[]
    for s in get_systems(node)
        isconnector(s) || continue
        push!(connectors, nameof(s))
    end
    if length(connectors) == 1
        return only(connectors)
    elseif isconnector(node)
        return nothing
    else
        throw(
            ArgumentError(
            lazy"""Cannot infer ports for node `$(nameof(node))`: it has $(length(connectors)) connector subsystems. Provide a `portmap` to `multiconnect` to select which port each edge uses. Available connectors: $(join(connectors, ", "))."""
        )
        )
    end
end

function _normalize_port_pair(pair, i, j)
    pair isa Union{Pair, Tuple{Any, Any}} ||
        throw(
            ArgumentError(
            lazy"Port map for edge $i -> $j must be a `Pair` or 2-`Tuple` of port names. Got `$pair`."
        )
        )
    src_port, dst_port = pair
    (src_port === nothing || src_port isa Symbol) ||
        throw(ArgumentError(lazy"Invalid source port `$src_port` for edge $i -> $j."))
    (dst_port === nothing || dst_port isa Symbol) ||
        throw(ArgumentError(lazy"Invalid destination port `$dst_port` for edge $i -> $j."))
    return (src_port, dst_port)
end

"""
    $(TYPEDSIGNATURES)

Resolve the `(src_port, dst_port)` pair for edge `i -> j` of a network using
`portmap`. `portmap` may be `nothing` (use each node's only connector), a `Pair` or
2-`Tuple` of port names applied to every edge, an `AbstractDict` keyed by
`Graphs.Edge` or `Tuple{Int, Int}` (both orientations are tried), or a callable
`portmap(i, j) -> Pair`.
"""
function _resolve_portmap(portmap::Union{Pair, Tuple{Any, Any}}, nodes, i::Int, j::Int)
    return _normalize_port_pair((portmap[1], portmap[2]), i, j)
end

function _resolve_portmap(portmap::AbstractDict, nodes, i::Int, j::Int)
    for key in (Graphs.Edge(i, j), Graphs.Edge(j, i), (i, j), (j, i))
        haskey(portmap, key) && return _normalize_port_pair(portmap[key], i, j)
    end
    throw(ArgumentError(lazy"No port mapping provided for edge $i -> $j."))
end

function _resolve_portmap(::Nothing, nodes, i::Int, j::Int)
    return (_default_network_port(nodes[i]), _default_network_port(nodes[j]))
end

function _resolve_portmap(portmap, nodes, i::Int, j::Int)
    return _normalize_port_pair(portmap(i, j), i, j)
end

"""
    multiconnect(nodes::AbstractVector{<:AbstractSystem}, edges::AbstractVector{ConnectionEdge})
    multiconnect(nodes::AbstractVector{<:AbstractSystem}, g::Graphs.AbstractGraph, portmap = nothing)

Create a single `Connection` equation connecting an array of node subsystems
according to a graph of labeled edges.

`nodes` is a vector of instantiated subsystems; different nodes may be
instantiations of different component forms, as long as the ports referenced by
each edge exist on the corresponding nodes. `edges` is a list of
[`ConnectionEdge`](@ref)s specifying, for each edge, which port of the source node
connects to which port of the destination node. Alternatively, a
`Graphs.AbstractGraph` can be supplied together with `portmap`, which selects the
port names used on each end of each edge. `portmap` may be:

- a `Pair` or 2-`Tuple` of port names, used for every edge (e.g. `:n => :p`),
- an `AbstractDict` mapping `Graphs.Edge` or `(i, j)` keys to such pairs,
- a callable `portmap(i, j) -> (src_port, dst_port)`, or
- `nothing`, in which case each node must have exactly one connector subsystem
  (or be a connector itself).

Each edge is lowered exactly like a pairwise [`connect`](@ref) between the two
ports. Acausal connector ports generate undirected (equality + sum-to-zero flow)
connection sets; causal input/output variable ports generate directed
`output ~ input` equations. Ports shared by multiple edges merge into a single
connection set, exactly as if the corresponding `connect` statements had been
written out.

Because the entire network is a single symbolic object, a system containing a
`multiconnect` equation stores the connection structure in O(1) additional
equations regardless of the number of nodes or edges; the materialized connection
equations are only generated during lowering.

# Example

```julia
# a chain of `n` two-port components, `comp[i].n` connected to `comp[i+1].p`
g = Graphs.path_graph(n)
eqs = [multiconnect(comps, g, :n => :p)]
```
"""
function multiconnect(
        nodes::AbstractVector{<:AbstractSystem}, edges::AbstractVector{ConnectionEdge}
)
    isempty(nodes) &&
        throw(ArgumentError("`multiconnect` requires at least one node."))
    _nodes = System[node for node in nodes]
    allunique(nameof, _nodes) ||
        throw(ArgumentError("`multiconnect` requires all nodes to have unique names."))
    n = length(_nodes)
    for edge in edges
        for (idx, side) in ((edge.src, "source"), (edge.dst, "destination"))
            1 <= idx <= n || throw(
                ArgumentError(
                lazy"Edge $edge references $side node $idx, but only $n nodes were given."
            )
            )
        end
        edge.src == edge.dst && edge.src_port == edge.dst_port &&
            throw(ArgumentError(lazy"Self-edge $edge connects a port to itself."))
        src_port = _network_port(_nodes[edge.src], edge.src_port)
        dst_port = _network_port(_nodes[edge.dst], edge.dst_port)
        (src_port isa AbstractSystem) == (dst_port isa AbstractSystem) || throw(
            ArgumentError(
            lazy"Edge $edge connects ports of different kinds: `$src_port` and `$dst_port`."
        )
        )
        if src_port isa SymbolicT
            validate_causal_variables_connection(SymbolicT[src_port, dst_port])
        end
    end
    network = ConnectionNetwork(_nodes, Vector{ConnectionEdge}(edges))
    return Equation(Connection(), Connection(network))
end

function multiconnect(
        nodes::AbstractVector{<:AbstractSystem}, g::Graphs.AbstractGraph, portmap = nothing
)
    Graphs.nv(g) == length(nodes) || throw(
        ArgumentError(
        lazy"Graph has $(Graphs.nv(g)) vertices, but $(length(nodes)) nodes were given."
    )
    )
    edges = ConnectionEdge[]
    sizehint!(edges, Graphs.ne(g))
    for e in Graphs.edges(g)
        i, j = Graphs.src(e), Graphs.dst(e)
        src_port, dst_port = _resolve_portmap(portmap, nodes, i, j)
        push!(edges, ConnectionEdge(i, j, src_port, dst_port))
    end
    return multiconnect(nodes, edges)
end

"""
    $(TYPEDSIGNATURES)

Generate connection sets for a [`ConnectionNetwork`](@ref). Each edge expands to
the same hyperedges a pairwise `connect` between the two resolved ports would
produce, so multi-edge ports merge into single connection sets through the
ordinary machinery.
"""
function _generate_connectionsets!(
        connection_state::AbstractConnectionState,
        namespace::Vector{Symbol}, network::ConnectionNetwork, isouter::IsOuter
)
    for edge in network.edges
        src_port = _network_port(network.nodes[edge.src], edge.src_port)
        dst_port = _network_port(network.nodes[edge.dst], edge.dst_port)
        if src_port isa AbstractSystem && dst_port isa AbstractSystem
            _generate_connectionsets!(
                connection_state, namespace, System[src_port, dst_port], isouter
            )
        else
            _generate_connectionsets!(
                connection_state, namespace, SymbolicT[src_port, dst_port], isouter
            )
        end
    end
    return nothing
end

"""
    $(TYPEDSIGNATURES)

Convert `x` to a tuple of unwrapped symbolic port variables.
"""
_port_var_tuple(x::Union{Tuple, AbstractVector}) = Tuple(unwrap(v) for v in x)
_port_var_tuple(x) = (unwrap(x),)

"""
    $(TYPEDSIGNATURES)

Validate that `var` is a one-dimensional symbolic array variable of length `n`
(the number of nodes in the network). Return the unwrapped variable.
"""
function _network_port_var(var, n::Int)
    var = unwrap(var)
    var isa SymbolicT ||
        throw(ArgumentError(lazy"Port `$var` is not a symbolic variable."))
    sh = SU.shape(var)
    sh isa SU.ShapeVecT && length(sh) == 1 || throw(
        ArgumentError(
        lazy"Port variable `$var` must be a one-dimensional symbolic array indexed by node."
    )
    )
    length(only(sh)) == n || throw(
        ArgumentError(
        lazy"Port variable `$var` has length $(length(only(sh))), expected `nv(graph) = $n`."
    )
    )
    return var::SymbolicT
end

"""
    $(TYPEDSIGNATURES)

Validate a causal port-variable pair for a network. Exactly one of `a`, `b` must
be an output and the other an input; `srcs`/`dsts` are the edge endpoints where
`a` is indexed by `srcs` and `b` by `dsts`. Each input port may be driven by at
most one edge.
"""
function _validate_causal_network_vars(a::SymbolicT, b::SymbolicT, srcs, dsts)
    driven = if isoutput(a) && !isinput(a) && isinput(b) && !isoutput(b)
        dsts
    elseif isinput(a) && !isoutput(a) && isoutput(b) && !isinput(b)
        srcs
    else
        throw(
            ArgumentError(
            lazy"Causal `multiconnect` requires each edge to pair an output variable with an input variable. Got `$a` (input = $(isinput(a)), output = $(isoutput(a))) and `$b` (input = $(isinput(b)), output = $(isoutput(b)))."
        )
        )
    end
    seen = Int[]
    for node in driven
        node in seen && throw(
            ArgumentError(
            lazy"Input port variable at node $node is driven by multiple edges in the network."
        )
        )
        push!(seen, node)
    end
    return nothing
end

"""
    $(TYPEDSIGNATURES)

Generate conservation equations for a flow port-variable pair `(a, b)` over the
edges `srcs`/`dsts` of a graph with `n` nodes. Each connected component ("net") of
the port graph produces one sum-to-zero equation over its member flow variables;
ports not touched by any edge produce singleton `f ~ 0` equations, matching the
semantics of unconnected acausal connectors. If `a` and `b` are the same variable
the port of a node is shared by all of its edges; otherwise the two sides are
distinct ports.
"""
function _flow_net_equations(a::SymbolicT, b::SymbolicT, n::Int, srcs, dsts)
    shared = isequal(a, b)
    nverts = shared ? n : 2n
    # union-find over port vertices
    parents = collect(1:nverts)
    function findroot(x::Int)
        while parents[x] != x
            parents[x] = parents[parents[x]]
            x = parents[x]
        end
        return x
    end
    for (i, j) in zip(srcs, dsts)
        ra = findroot(i)
        rb = findroot(shared ? j : n + j)
        ra != rb && (parents[rb] = ra)
    end
    nets = Dict{Int, Vector{SymbolicT}}()
    for v in 1:nverts
        var = v <= n ? a[v] : b[v - n]
        members = get!(nets, findroot(v), SymbolicT[])
        push!(members, var)
    end
    eqs = Equation[]
    sizehint!(eqs, length(nets))
    for root in sort!(collect(keys(nets)))
        members = nets[root]
        push!(eqs, Symbolics.COMMON_ZERO ~ SU.add_worker(VartypeT, members))
    end
    return eqs
end

"""
    multiconnect(portspec::Pair, g::Graphs.AbstractGraph)

Compact "array of components" form of [`multiconnect`](@ref), where every node is
an index into shared symbolic array variables rather than a distinct subsystem.

Each side of `portspec` is a symbolic array variable (or a tuple/vector of such
variables) of length `nv(g)`, giving the per-node port variables on the source
and destination side of each edge. For every edge `i -> j`, the k-th source port
variable at node `i` is connected to the k-th destination port variable at node
`j`.

Returns a `Vector{Equation}` of *compact* array equations - the connection graph
is encoded in sparse selection matrices, so the number of generated symbolic
equations is independent of the number of edges:

- `Equality` (default) port variables `a`, `b` generate `Ps * a ~ Pd * b`, where
  `Ps`/`Pd` are sparse matrices selecting the source/destination endpoint of each
  edge. This represents all `a[i] ~ b[j]` edge equalities in a single symbolic
  equation.
- `Flow` port variables generate one sum-to-zero conservation equation per
  connected component (net) of the port graph, plus `f[i] ~ 0` for ports not
  touched by any edge - the same equations acausal connection expansion produces.
- Variables with `[input = true]`/`[output = true]` metadata are treated as
  causal: each edge must pair an output with an input, and each input may be
  driven by at most one edge.

!!! note
    Unlike the subsystem form, all edges referencing a given flow variable must
    be supplied in a single `multiconnect` call for the net conservation
    equations to be correct.

# Example

```julia
@variables v(t)[1:n] i(t)[1:n]  # per-node potential and flow
eqs = multiconnect((v, i) => (v, i), g)  # Kirchhoff network on graph `g`
```
"""
function multiconnect(portspec::Pair, g::Graphs.AbstractGraph)
    avars = _port_var_tuple(portspec.first)
    bvars = _port_var_tuple(portspec.second)
    length(avars) == length(bvars) || throw(
        ArgumentError(
        lazy"Source and destination port lists must have the same length. Got $(length(avars)) and $(length(bvars))."
    )
    )
    isempty(avars) &&
        throw(ArgumentError("`multiconnect` requires at least one port variable."))
    n = Graphs.nv(g)
    n == 0 && throw(ArgumentError("`multiconnect` requires a non-empty graph."))
    srcs = Int[]
    dsts = Int[]
    sizehint!(srcs, Graphs.ne(g))
    sizehint!(dsts, Graphs.ne(g))
    for e in Graphs.edges(g)
        push!(srcs, Graphs.src(e))
        push!(dsts, Graphs.dst(e))
    end
    nedges = length(srcs)
    Ps = sparse(1:nedges, srcs, 1.0, nedges, n)
    Pd = sparse(1:nedges, dsts, 1.0, nedges, n)
    eqs = Equation[]
    for (a, b) in zip(avars, bvars)
        a = _network_port_var(a, n)
        b = _network_port_var(b, n)
        if (isinput(a) || isoutput(a)) || (isinput(b) || isoutput(b))
            _validate_causal_network_vars(a, b, srcs, dsts)
            push!(eqs, Ps * a ~ Pd * b)
            continue
        end
        ctype_a = get_connection_type(a)
        ctype_b = get_connection_type(b)
        ctype_a === ctype_b || throw(
            ArgumentError(
            lazy"Connected port variables `$a` and `$b` have different connect types `$ctype_a` and `$ctype_b`."
        )
        )
        if ctype_a === Flow
            append!(eqs, _flow_net_equations(a, b, n, srcs, dsts))
        elseif ctype_a === Stream
            throw(
                ArgumentError(
                "`Stream` port variables are not supported by `multiconnect`."
            )
            )
        else # Equality
            push!(eqs, Ps * a ~ Pd * b)
        end
    end
    return eqs
end
