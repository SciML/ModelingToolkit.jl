"""
    $(TYPEDEF)

Struct used to represent a connection equation. A connection equation is an `Equation`
where the LHS is an empty `Connection(nothing)` and the RHS is a `Connection` containing
the connected connectors.

For special types of connections, the LHS `Connection` can contain relevant metadata.
"""
struct Connection
    systems::Any
end

Base.broadcastable(x::Connection) = Ref(x)
Connection() = Connection(nothing)
Base.hash(c::Connection, seed::UInt) = hash(c.systems, (0xc80093537bdc1311 % UInt) ⊻ seed)
Symbolics.hide_lhs(_::Connection) = true

"""
    $(TYPEDSIGNATURES)

Connect multiple connectors created via `@connector`. All connected connectors
must be unique.
"""
function connect(sys1::AbstractSystem, sys2::AbstractSystem, syss::AbstractSystem...)
    syss = (sys1, sys2, syss...)
    length(unique(nameof, syss)) == length(syss) || error("connect takes distinct systems!")
    Equation(Connection(), Connection(syss)) # the RHS are connected systems
end

const _debug_mode = Base.JLOptions().check_bounds == 1

function Base.show(io::IO, c::Connection)
    print(io, "connect(")
    if c.systems isa AbstractArray || c.systems isa Tuple
        n = length(c.systems)
        for (i, s) in enumerate(c.systems)
            str = join(split(string(nameof(s)), NAMESPACE_SEPARATOR), '.')
            print(io, str)
            i != n && print(io, ", ")
        end
    end
    print(io, ")")
end

@latexrecipe function f(c::Connection)
    index --> :subscript
    return Expr(:call, :connect, map(nameof, c.systems)...)
end

function Base.show(io::IO, ::MIME"text/latex", ap::Connection)
    print(io, latexify(ap))
end

isconnection(_) = false
isconnection(_::Connection) = true

"""
    $(TYPEDSIGNATURES)

Adds a domain only connection equation, through and across state equations are not generated.
"""
function domain_connect(sys1::AbstractSystem, sys2::AbstractSystem, syss::AbstractSystem...)
    syss = (sys1, sys2, syss...)
    length(unique(nameof, syss)) == length(syss) || error("connect takes distinct systems!")
    Equation(Connection(:domain), Connection(syss)) # the RHS are connected systems
end

"""
    $(TYPEDSIGNATURES)

Get the connection type of symbolic variable `s` from the `VariableConnectType` metadata.
Defaults to `Equality` if not present.
"""
function get_connection_type(s::Symbolic)
    s = unwrap(s)
    if iscall(s) && operation(s) === getindex
        s = arguments(s)[1]
    end
    getmetadata(s, VariableConnectType, Equality)
end

"""
    $(TYPEDSIGNATURES)

Mark a system constructor function as building a connector. For example,

```julia
@connector function ElectricalPin(; name, v = nothing, i = nothing)
    @variables begin
        v(t) = v, [description = "Potential at the pin [V]"]
        i(t) = i, [connect = Flow, description = "Current flowing into the pin [A]"]
    end
    return System(Equation[], t, [v, i], []; name)
end
```

Since connectors only declare variables, the equivalent shorthand syntax can also be used:

```julia
@connector Pin begin
    v(t), [description = "Potential at the pin [V]"]
    i(t), [connect = Flow, description = "Current flowing into the pin [A]"]
end
```

ModelingToolkit systems are either components or connectors. Components define dynamics of
the model. Connectors are used to connect components together. See the
[Model building reference](@ref model_building_api) section of the documentation for more
information.

See also: [`@component`](@ref).
"""
macro connector(expr)
    esc(component_post_processing(expr, true))
end

abstract type AbstractConnectorType end
struct StreamConnector <: AbstractConnectorType end
struct RegularConnector <: AbstractConnectorType end
struct DomainConnector <: AbstractConnectorType end

"""
    $(TYPEDSIGNATURES)

Return an `AbstractConnectorType` denoting the type of connector that `sys` is.
Domain connectors have a single `Flow` unknown. Stream connectors have a single
`Flow` variable and multiple `Stream` variables. Any other type of connector is
a "regular" connector.
"""
function connector_type(sys::AbstractSystem)
    unkvars = get_unknowns(sys)
    n_stream = 0
    n_flow = 0
    n_regular = 0 # unknown that is not input, output, stream, or flow.
    for s in unkvars
        vtype = get_connection_type(s)
        if vtype === Stream
            isarray(s) && error("Array stream variables are not supported. Got $s.")
            n_stream += 1
        elseif vtype === Flow
            n_flow += 1
        elseif !(isinput(s) || isoutput(s))
            n_regular += 1
        end
    end
    (n_stream > 0 && n_flow > 1) &&
        error("There are multiple flow variables in the stream connector $(nameof(sys))!")
    if n_flow == 1 && length(unkvars) == 1
        return DomainConnector()
    end
    if n_flow != n_regular && !isframe(sys)
        @warn "$(nameof(sys)) contains $n_flow flow variables, yet $n_regular regular " *
              "(non-flow, non-stream, non-input, non-output) variables. " *
              "This could lead to imbalanced model that are difficult to debug. " *
              "Consider marking some of the regular variables as input/output variables."
    end
    n_stream > 0 ? StreamConnector() : RegularConnector()
end

is_domain_connector(s) = isconnector(s) && get_connector_type(s) === DomainConnector()

get_systems(c::Connection) = c.systems

"""
    $(TYPEDSIGNATURES)

`instream` is used when modeling stream connections. It is only allowed to be used on
`Stream` variables.

Refer to the [Connection semantics](@ref connect_semantics) section of the docs for more
information.
"""
instream(a) = term(instream, unwrap(a), type = symtype(a))
SymbolicUtils.promote_symtype(::typeof(instream), _) = Real

isconnector(s::AbstractSystem) = has_connector_type(s) && get_connector_type(s) !== nothing

"""
    $(TYPEDEF)

Utility struct which wraps a symbolic variable used in a `Connection` to enable `Base.show`
to work.
"""
struct SymbolicWithNameof
    var::Any
end

function Base.nameof(x::SymbolicWithNameof)
    return Symbol(x.var)
end

is_causal_variable_connection(c) = false
function is_causal_variable_connection(c::Connection)
    all(x -> x isa SymbolicWithNameof, get_systems(c))
end

const ConnectableSymbolicT = Union{BasicSymbolic, Num, Symbolics.Arr}

function NonCausalVariableError(vars)
    names = join(map(var -> "  " * string(var), vars), "\n")
    ArgumentError("""
    Only causal variables can be used in a `connect` statement. Each variable must be \
    either an input or an output. Mark a variable as input using the `[input = true]` \
    variable metadata or as an output using the `[output = true]` variable metadata.

    The following variables were found to be non-causal:
    $names
    """)
end

"""
    $(TYPEDSIGNATURES)

Perform validation for a connect statement involving causal variables.
"""
function validate_causal_variables_connection(allvars)
    var1 = allvars[1]
    var2 = allvars[2]
    vars = Base.tail(Base.tail(allvars))
    for var in allvars
        vtype = getvariabletype(var)
        vtype === VARIABLE ||
            throw(ArgumentError("Expected $var to be of kind `$VARIABLE`. Got `$vtype`."))
    end
    if length(unique(allvars)) !== length(allvars)
        throw(ArgumentError("Expected all connection variables to be unique. Got variables $allvars which contains duplicate entries."))
    end
    allsizes = map(size, allvars)
    if !allequal(allsizes)
        throw(ArgumentError("Expected all connection variables to have the same size. Got variables $allvars with sizes $allsizes respectively."))
    end
    non_causal_variables = filter(allvars) do var
        !isinput(var) && !isoutput(var)
    end
    isempty(non_causal_variables) || throw(NonCausalVariableError(non_causal_variables))
end

"""
    $(TYPEDSIGNATURES)

Connect multiple causal variables. The first variable must be an output, and all subsequent
variables must be inputs. The statement `connect(var1, var2, var3, ...)` expands to:

```julia
var1 ~ var2
var1 ~ var3
# ...
```
"""
function connect(var1::ConnectableSymbolicT, var2::ConnectableSymbolicT,
        vars::ConnectableSymbolicT...)
    allvars = (var1, var2, vars...)
    validate_causal_variables_connection(allvars)
    return Equation(Connection(), Connection(map(SymbolicWithNameof, unwrap.(allvars))))
end

"""
    $(METHODLIST)

Add all `instream(..)` expressions to `set`.
"""
function collect_instream!(set, eq::Equation)
    collect_instream!(set, eq.lhs) | collect_instream!(set, eq.rhs)
end

function collect_instream!(set, expr, occurs = false)
    iscall(expr) || return occurs
    op = operation(expr)
    op === instream && (push!(set, expr); occurs = true)
    for a in SymbolicUtils.arguments(expr)
        occurs |= collect_instream!(set, a, occurs)
    end
    return occurs
end

#positivemax(m, ::Any; tol=nothing)= max(m, something(tol, 1e-8))
#_positivemax(m, tol) = ifelse((-tol <= m) & (m <= tol), ((3 * tol - m) * (tol + m)^3)/(16 * tol^3) + tol, max(m, tol))
function _positivemax(m, si)
    T = typeof(m)
    relativeTolerance = 1e-4
    nominal = one(T)
    eps = relativeTolerance * nominal
    alpha = if si > eps
        one(T)
    else
        if si > 0
            (si / eps)^2 * (3 - 2 * si / eps)
        else
            zero(T)
        end
    end
    alpha * max(m, 0) + (1 - alpha) * eps
end
@register_symbolic _positivemax(m, tol)
positivemax(m, ::Any; tol = nothing) = _positivemax(m, tol)
mydiv(num, den) =
    if den == 0
        error()
    else
        num / den
    end
@register_symbolic mydiv(n, d)

"""
    $(TYPEDSIGNATURES)

Return a function which checks whether the connector (system) passed to it is an outside
connector of `sys`. The function can also be given the name of a system as a `Symbol`.
"""
function generate_isouter(sys::AbstractSystem)
    outer_connectors = Symbol[]
    for s in get_systems(sys)
        n = nameof(s)
        isconnector(s) && push!(outer_connectors, n)
    end
    let outer_connectors = outer_connectors
        function isouter(sys)::Bool
            s = string(nameof(sys))
            isconnector(sys) || error("$s is not a connector!")
            idx = findfirst(isequal(NAMESPACE_SEPARATOR), s)
            parent_name = Symbol(idx === nothing ? s : s[1:prevind(s, idx)])
            isouter(parent_name)
        end
        function isouter(name::Symbol)::Bool
            return name in outer_connectors
        end
    end
end

@noinline function connection_error(ss)
    error("Different types of connectors are in one connection statement: <$(map(nameof, ss))>")
end

abstract type IsFrame end

"Return true if the system is a 3D multibody frame, otherwise return false."
function isframe(sys)
    getmetadata(sys, IsFrame, false)
end

abstract type FrameOrientation end

"Return orientation object of a multibody frame."
function ori(sys)
    getmetadata(sys, FrameOrientation, nothing)
end

"""
Connection type used in `ConnectionVertex` for a causal input variable. `I` is an object
that can be passed to `getindex` as an index denoting the index in the variable for
causal array variables. For non-array variables this should be `1`.
"""
abstract type InputVar{I} end
"""
Connection type used in `ConnectionVertex` for a causal output variable. `I` is an object
that can be passed to `getindex` as an index denoting the index in the variable for
causal array variables. For non-array variables this should be `1`.
"""
abstract type OutputVar{I} end

"""
    $(METHODLIST)

Get the contained index in an `InputVar` or `OutputVar` type.
"""
index_from_type(::Type{InputVar{I}}) where {I} = I
index_from_type(::Type{OutputVar{I}}) where {I} = I

"""
    $(TYPEDSIGNATURES)

Chain `getproperty` calls on sys in the order given by `names` and return the unwrapped
result.
"""
function iterative_getproperty(sys::AbstractSystem, names::AbstractVector{Symbol})
    # we don't want to namespace the first time
    result = toggle_namespacing(sys, false)
    for name in names
        result = getproperty(result, name)
    end
    return unwrap(result)
end

"""
    $(TYPEDSIGNATURES)

Return the variable/subsystem of `sys` referred to by vertex `vert`.
"""
function variable_from_vertex(sys::AbstractSystem, vert::ConnectionVertex)
    value = iterative_getproperty(sys, vert.name)
    value isa AbstractSystem && return value
    vert.type <: Union{InputVar, OutputVar} || return value
    # index possibly array causal variable
    unwrap(wrap(value)[index_from_type(vert.type)])
end

"""
    $(TYPEDSIGNATURES)

Given `connected`, the list of connected variables/systems, generate the appropriate
connection sets and add them to `connection_state`. Update both the connection network and
domain network as necessary. `namespace` is the path from the root system to the system in
which the [`connect`](@ref) equation containing `connected` is located. `isouter` is the
function returned from [`generate_isouter`](@ref) for the system referred to by
`namespace`.

`namespace` must not contain the name of the root system.
"""
function generate_connectionsets!(connection_state::AbstractConnectionState,
        namespace::Vector{Symbol}, connected, isouter)
    initial_len = length(namespace)
    _generate_connectionsets!(connection_state, namespace, connected, isouter)
    # Enforce postcondition as a sanity check that the namespacing is implemented correctly
    length(namespace) == initial_len || throw(NotPossibleError())
    return nothing
end

function _generate_connectionsets!(connection_state::AbstractConnectionState,
        namespace::Vector{Symbol},
        connected_vars::Union{
            AbstractVector{SymbolicWithNameof}, Tuple{Vararg{SymbolicWithNameof}}},
        isouter)
    # unwrap the `SymbolicWithNameof` into the contained symbolic variables.
    connected_vars = map(x -> x.var, connected_vars)
    _generate_connectionsets!(connection_state, namespace, connected_vars, isouter)
end

function _generate_connectionsets!(connection_state::AbstractConnectionState,
        namespace::Vector{Symbol},
        connected_vars::Union{
            AbstractVector{<:BasicSymbolic}, Tuple{Vararg{BasicSymbolic}}},
        isouter)
    # NOTE: variable connections don't populate the domain network

    # wrap to be able to call `eachindex` on a non-array variable
    representative = wrap(first(connected_vars))
    # all of them have the same size, but may have different axes/shape
    # so we iterate over `eachindex(eachindex(..))` since that is identical for all
    for sz_i in eachindex(eachindex(representative))
        hyperedge = map(connected_vars) do var
            var = unwrap(var)
            var_ns = namespace_hierarchy(getname(var))
            i = eachindex(wrap(var))[sz_i]

            is_input = isinput(var)
            is_output = isoutput(var)
            if is_input && is_output
                names = join(string.(connected_vars), ", ")
                throw(ArgumentError("""
                Variable $var in connection `connect($names)` is both input and output.
                """))
            elseif is_input
                type = InputVar{i}
            elseif is_output
                type = OutputVar{i}
            else
                names = join(string.(connected_vars), ", ")
                throw(ArgumentError("""
                Variable $var in connection `connect($names)` is neither input nor output.
                """))
            end

            return ConnectionVertex(
                [namespace; var_ns], length(var_ns) == 1 || isouter(var_ns[1]), type)
        end
        add_connection_edge!(connection_state, hyperedge)

        # Removed analysis points generate causal connections in the negative graph. These
        # should also remove `Equality` connections involving the same variables, so also
        # add an `Equality` variant of the edge.
        if connection_state isa NegativeConnectionState
            hyperedge = map(hyperedge) do cvert
                ConnectionVertex(cvert.name, cvert.isouter, Equality)
            end
            add_connection_edge!(connection_state, hyperedge)
        end
    end
end

function _generate_connectionsets!(connection_state::AbstractConnectionState,
        namespace::Vector{Symbol},
        systems::Union{AbstractVector{<:AbstractSystem}, Tuple{Vararg{AbstractSystem}}},
        isouter)
    regular_systems = System[]
    domain_system = nothing
    for s in systems
        if is_domain_connector(s)
            if domain_system === nothing
                domain_system = s
            else
                names = join(map(string ∘ nameof, systems), ",")
                error("connect($names) contains multiple source domain connectors. There can only be one!")
            end
        else
            push!(regular_systems, s)
        end
    end

    @assert !isempty(regular_systems)

    systems = regular_systems
    # There is a domain being connected here. In such a case, we only connect the
    # flow variable common between the domain setter and all other connectors in the
    # normal connection graph. The domain graph connects all these subsystems.
    if domain_system !== nothing
        hyperedge = ConnectionVertex[]
        domain_hyperedge = ConnectionVertex[]
        sizehint!(hyperedge, length(systems) + 1)
        sizehint!(domain_hyperedge, length(systems) + 1)

        dv = only(unknowns(domain_system))
        push!(namespace, nameof(domain_system))
        dv_vertex = ConnectionVertex(namespace, dv, false)
        domain_vertex = ConnectionVertex(namespace)
        pop!(namespace)

        push!(domain_hyperedge, domain_vertex)
        push!(hyperedge, dv_vertex)

        for (i, sys) in enumerate(systems)
            sts = unknowns(sys)
            sys_is_outer = isouter(sys)

            # add this system to the namespace so all vertices created from its unknowns
            # are properly namespaced
            sysname = nameof(sys)
            sys_ns = namespace_hierarchy(sysname)
            append!(namespace, sys_ns)
            for v in sts
                vtype = get_connection_type(v)
                # ignore all non-flow vertices in connectors
                (vtype === Flow && isequal(v, dv)) || continue

                vertex = ConnectionVertex(namespace, v, sys_is_outer)
                # vertices in the domain graph are systems with isouter=true and type=Flow
                sys_vertex = ConnectionVertex(namespace)
                push!(hyperedge, vertex)
                push!(domain_hyperedge, sys_vertex)
            end
            # remember to remove the added namespace!
            foreach(_ -> pop!(namespace), sys_ns)
        end
        @assert length(hyperedge) > 1
        @assert length(domain_hyperedge) == length(hyperedge)

        add_connection_edge!(connection_state, hyperedge)
        add_domain_connection_edge!(connection_state, domain_hyperedge)
        return
    end
    sys1 = first(systems)
    sys1_dvs = unknowns(sys1)
    # Add 9 orientation variables if connection is between multibody frames
    if isframe(sys1) # Multibody
        O = ori(sys1)
        orientation_vars = Symbolics.unwrap.(collect(vec(O.R)))
        sys1_dvs = [sys1_dvs; orientation_vars]
    end
    sys1_dvs_set = Set(sys1_dvs)
    num_unknowns = length(sys1_dvs)

    # We first build sets of all vertices that are connected together
    var_sets = [ConnectionVertex[] for _ in 1:num_unknowns]
    domain_hyperedge = ConnectionVertex[]
    for (i, sys) in enumerate(systems)
        unknown_vars = unknowns(sys)
        # Add 9 orientation variables if connection is between multibody frames
        if isframe(sys) # Multibody
            O = ori(sys)
            orientation_vars = Symbolics.unwrap.(vec(O.R))
            unknown_vars = [unknown_vars; orientation_vars]
        end
        # Error if any subsequent systems do not have the same number of unknowns
        # or have unknowns not in the others.
        if i != 1 &&
           (num_unknowns != length(unknown_vars) || any(!in(sys1_dvs_set), unknown_vars))
            connection_error(systems)
        end
        # add this system to the namespace so all vertices created from its unknowns
        # are properly namespaced
        sysname = nameof(sys)
        sys_ns = namespace_hierarchy(sysname)
        append!(namespace, sys_ns)
        sys_is_outer = isouter(sys)
        for (j, v) in enumerate(unknown_vars)
            push!(var_sets[j], ConnectionVertex(namespace, v, sys_is_outer))
        end
        domain_vertex = ConnectionVertex(namespace)
        push!(domain_hyperedge, domain_vertex)
        # remember to remove the added namespace!
        foreach(_ -> pop!(namespace), sys_ns)
    end
    for var_set in var_sets
        # all connected variables should have the same type
        if !allequal(Iterators.map(cvert -> cvert.type, var_set))
            connection_error(systems)
        end
        # add edges
        add_connection_edge!(connection_state, var_set)
    end
    add_domain_connection_edge!(connection_state, domain_hyperedge)
end

"""
    $(TYPEDSIGNATURES)

Generate the merged connection sets and connected domain sets for system `sys`. Also
removes all `connect` equations in `sys`. Return the modified system and a tuple of the
connection sets and domain sets. Also scalarizes array equations in the system.
"""
function generate_connection_set(sys::AbstractSystem)
    # generate the states
    connection_state = ConnectionState()
    negative_connection_state = NegativeConnectionState()
    # the root system isn't added to the namespace, which we handle by not namespacing it
    sys = toggle_namespacing(sys, false)
    sys = generate_connection_set!(
        connection_state, negative_connection_state, sys, Symbol[])
    remove_negative_connections!(connection_state, negative_connection_state)

    return sys, connectionsets(connection_state)
end

"""
    $(TYPEDSIGNATURES)

Appropriately handle the equation `eq` depending on whether it is a normal or connection
equation. For normal equations, it is expected that `eqs` is a buffer to which the equation
can be pushed, unmodified. Connection equations update the given `state`. The equation is
present at the path in the hierarchical system given by `namespace`. `isouter` is the
function returned from `generate_isouter`.
"""
function handle_maybe_connect_equation!(eqs, state::AbstractConnectionState,
        eq::Equation, namespace::Vector{Symbol}, isouter)
    lhs = eq.lhs
    rhs = eq.rhs

    if !(lhs isa Connection)
        # split connections and equations
        if eq.lhs isa AbstractArray || eq.rhs isa AbstractArray
            append!(eqs, Symbolics.scalarize(eq))
        else
            push!(eqs, eq)
        end
        return
    end

    if get_systems(lhs) === :domain
        # This is a domain connection, so we only update the domain connection graph
        hyperedge = map(get_systems(rhs)) do sys
            sys isa AbstractSystem || error("Domain connections can only connect systems!")
            sysname = nameof(sys)
            sys_ns = namespace_hierarchy(sysname)
            append!(namespace, sys_ns)
            vertex = ConnectionVertex(namespace)
            foreach(_ -> pop!(namespace), sys_ns)
            return vertex
        end
        add_domain_connection_edge!(state, hyperedge)
    else
        connected_systems = get_systems(rhs)
        generate_connectionsets!(state, namespace, connected_systems, isouter)
    end
    return nothing
end

"""
    $(TYPEDSIGNATURES)

Generate the appropriate connection sets from `connect` equations present in the
hierarchical system `sys`. This is a recursive function that descends the hierarchy. If
`sys` is the root system, then `does_namespacing(sys)` must be `false` and `namespace`
should be empty. It is essential that the traversal is preorder.

## Arguments

- `connection_state`: The connection state keeping track of the connection network and the
  domain network.
- `negative_connection_state`: The connection state that tracks connections removed by
  analysis point transformations. These removed connections are stored in the
  `ignored_connections` field of the system.
- `namespace`: The path of names from the root system to the current system. This should
  not include the name of the root system.
"""
function generate_connection_set!(connection_state::ConnectionState,
        negative_connection_state::NegativeConnectionState,
        sys::AbstractSystem, namespace::Vector{Symbol})
    initial_len = length(namespace)
    res = _generate_connection_set!(
        connection_state, negative_connection_state, sys, namespace)
    # Enforce postcondition as a sanity check that the recursion is implemented correctly
    length(namespace) == initial_len || throw(NotPossibleError())
    return res
end

function _generate_connection_set!(connection_state::ConnectionState,
        negative_connection_state::NegativeConnectionState,
        sys::AbstractSystem, namespace::Vector{Symbol})
    # This function recurses down the system tree. Each system adds its name and pops
    # it before returning. We don't add the root system, which is handled by assuming
    # it doesn't do namespacing.
    does_namespacing(sys) && push!(namespace, nameof(sys))
    subsys = get_systems(sys)

    isouter = generate_isouter(sys)
    eqs′ = get_eqs(sys)
    eqs = Equation[]

    # generate connection equations and separate out non-connection equations
    for eq in eqs′
        handle_maybe_connect_equation!(eqs, connection_state, eq, namespace, isouter)
    end

    # go through the removed connections and update the negative graph
    for conn in something(get_ignored_connections(sys), ())
        eq = Equation(Connection(), conn)
        # there won't be any standard equations, so we can pass `nothing` instead of
        # `eqs`.
        handle_maybe_connect_equation!(
            nothing, negative_connection_state, eq, namespace, isouter)
    end

    # all connectors are eventually inside connectors, and all flow variables
    # need at least a singleton connectionset (hyperedge) with the inside variant
    for s in subsys
        isconnector(s) || continue
        is_domain_connector(s) && continue
        push!(namespace, nameof(s))
        for v in unknowns(s)
            Flow === get_connection_type(v) || continue
            add_connection_edge!(connection_state, (ConnectionVertex(namespace, v, false),))
        end
        pop!(namespace)
    end

    # recurse down the hierarchy
    @set! sys.systems = map(subsys) do s
        generate_connection_set!(connection_state, negative_connection_state, s, namespace)
    end
    @set! sys.eqs = eqs
    # Remember to pop the name at the end!
    does_namespacing(sys) && pop!(namespace)
    return sys
end

"""
    $(TYPEDSIGNATURES)

Generate connection equations for the connection sets given by `csets`. This does not
handle stream connections. Return the generated equations and the stream connection sets.
"""
function generate_connection_equations_and_stream_connections(
        sys::AbstractSystem, csets::Vector{Vector{ConnectionVertex}})
    eqs = Equation[]
    stream_connections = Vector{ConnectionVertex}[]

    for cset in csets
        cvert = cset[1]
        var = variable_from_vertex(sys, cvert)::BasicSymbolic
        vtype = cvert.type
        if vtype <: Union{InputVar, OutputVar}
            length(cset) > 1 || continue
            inner_output = nothing
            outer_input = nothing
            for cvert in cset
                if cvert.isouter && cvert.type <: InputVar
                    if outer_input !== nothing
                        error("""
                        Found two outer input connectors `$outer_input` and `$cvert` in the
                        same connection set.
                        """)
                    end
                    outer_input = cvert
                elseif !cvert.isouter && cvert.type <: OutputVar
                    if inner_output !== nothing
                        error("""
                        Found two inner output connectors `$inner_output` and `$cvert` in
                        the same connection set.
                        """)
                    end
                    inner_output = cvert
                end
            end
            root_vert = something(inner_output, outer_input)
            root_var = variable_from_vertex(sys, root_vert)
            for cvert in cset
                isequal(cvert, root_vert) && continue
                push!(eqs, variable_from_vertex(sys, cvert) ~ root_var)
            end
        elseif vtype === Stream
            push!(stream_connections, cset)
        elseif vtype === Flow
            # arrays have to be broadcasted to be added/subtracted/negated which leads
            # to bad-looking equations. Just generate scalar equations instead since
            # mtkcompile will scalarize anyway.
            representative = variable_from_vertex(sys, cset[1])
            # each variable can have different axes, but they all have the same size
            for sz_i in eachindex(eachindex(wrap(representative)))
                rhs = 0
                for cvert in cset
                    # all of this wrapping/unwrapping is necessary because the relevant
                    # methods are defined on `Arr/Num` and not `BasicSymbolic`.
                    v = variable_from_vertex(sys, cvert)::BasicSymbolic
                    idxs = eachindex(wrap(v))
                    v = unwrap(wrap(v)[idxs[sz_i]])
                    rhs += cvert.isouter ? unwrap(-wrap(v)) : v
                end
                push!(eqs, 0 ~ rhs)
            end
        else # Equality
            vars = map(Base.Fix1(variable_from_vertex, sys), cset)
            outer_input = inner_output = nothing
            all_io = true
            # attempt to interpret the equality as a causal connectionset if
            # possible
            for (cvert, vert) in zip(cset, vars)
                is_i = isinput(vert)
                is_o = isoutput(vert)
                all_io &= is_i || is_o
                all_io || break
                if cvert.isouter && is_i && outer_input === nothing
                    outer_input = cvert
                elseif !cvert.isouter && is_o && inner_output === nothing
                    inner_output = cvert
                end
            end
            # this doesn't necessarily mean this is a well-structured causal connection,
            # but it is sufficient and we're generating equalities anyway.
            if all_io && xor(outer_input !== nothing, inner_output !== nothing)
                root_vert = something(inner_output, outer_input)
                root_var = variable_from_vertex(sys, root_vert)
                for (cvert, var) in zip(cset, vars)
                    isequal(cvert, root_vert) && continue
                    push!(eqs, var ~ root_var)
                end
            else
                base = variable_from_vertex(sys, cset[1])
                for i in 2:length(cset)
                    v = vars[i]
                    push!(eqs, base ~ v)
                end
            end
        end
    end
    eqs, stream_connections
end

"""
    $(TYPEDSIGNATURES)

Generate the defaults for parameters in the domain sets given by `domain_csets`.
"""
function domain_defaults(
        sys::AbstractSystem, domain_csets::Vector{Vector{ConnectionVertex}})
    defs = Dict()
    for cset in domain_csets
        systems = map(Base.Fix1(variable_from_vertex, sys), cset)
        @assert all(x -> x isa AbstractSystem, systems)
        idx = findfirst(is_domain_connector, systems)
        idx === nothing && continue
        domain_sys = systems[idx]
        # note that these will not be namespaced with `domain_sys`.
        domain_defs = defaults(domain_sys)
        for (j, csys) in enumerate(systems)
            j == idx && continue
            if is_domain_connector(csys)
                throw(ArgumentError("""
                Domain sources $(nameof(domain_sys)) and $(nameof(csys)) are connected!
                """))
            end
            for par in parameters(csys)
                defval = get(domain_defs, par, nothing)
                defval === nothing && continue
                defs[parameters(csys, par)] = parameters(domain_sys, par)
            end
        end
    end
    return defs
end

"""
    $(TYPEDSIGNATURES)

Given a hierarchical system with [`connect`](@ref) equations, expand the connection
equations and return the new system. `tol` is the tolerance for handling the singularities
in stream connection equations that happen when a flow variable approaches zero.
"""
function expand_connections(sys::AbstractSystem; tol = 1e-10)
    # turn analysis points into standard connection equations
    sys = remove_analysis_points(sys)
    # generate the connection sets
    sys, (csets, domain_csets) = generate_connection_set(sys)
    # generate equations, and stream equations
    ceqs, instream_csets = generate_connection_equations_and_stream_connections(sys, csets)
    stream_eqs, instream_subs = expand_instream(instream_csets, sys; tol = tol)

    eqs = [equations(sys); ceqs; stream_eqs]
    if !isempty(instream_subs)
        # substitute `instream(..)` expressions with their new values
        for i in eachindex(eqs)
            eqs[i] = fixpoint_sub(
                eqs[i], instream_subs; maxiters = max(length(instream_subs), 10))
        end
    end
    # get the defaults for domain networks
    d_defs = domain_defaults(sys, domain_csets)
    # build the new system
    sys = flatten(sys, true)
    @set! sys.eqs = eqs
    @set! sys.defaults = merge(get_defaults(sys), d_defs)
end

"""
    $(TYPEDSIGNATURES)

Given a connection vertex `cvert` referring to a variable in a connector in `sys`, return
the flow variable in that connector.
"""
function get_flowvar(sys::AbstractSystem, cvert::ConnectionVertex)
    parent_names = @view cvert.name[1:(end - 1)]
    parent_sys = iterative_getproperty(sys, parent_names)
    for var in unknowns(parent_sys)
        type = get_connection_type(var)
        type == Flow || continue
        return unwrap(unknowns(parent_sys, var))
    end
    throw(ArgumentError("There is no flow variable in system `$(nameof(parent_sys))`"))
end

"""
    $(TYPEDSIGNATURES)

Given connection sets of stream variables in `sys`, return the additional equations to add
to the system and the substitutions to make to handle `instream(..)` expressions. `tol` is
the tolerance for handling singularities in stream connection equations when the flow
variable approaches zero.
"""
function expand_instream(csets::Vector{Vector{ConnectionVertex}}, sys::AbstractSystem;
        tol = 1e-8)
    eqs = equations(sys)
    # collect all `instream` terms in the equations
    instream_exprs = Set{BasicSymbolic}()
    for eq in eqs
        collect_instream!(instream_exprs, eq)
    end

    # specifically substitute `instream(x[i]) => instream(x)[i]`
    instream_subs = Dict{BasicSymbolic, BasicSymbolic}()
    for expr in instream_exprs
        stream_var = only(arguments(expr))
        iscall(stream_var) && operation(stream_var) === getindex || continue
        args = arguments(stream_var)
        new_expr = Symbolics.array_term(
            instream, args[1]; size = size(args[1]), ndims = ndims(args[1]))[args[2:end]...]
        instream_subs[expr] = new_expr
    end

    # for all the newly added `instream(x)[i]`, add `instream(x)` to `instream_exprs`
    # also remove all `instream(x[i])`
    for (k, v) in instream_subs
        push!(instream_exprs, arguments(v)[1])
        delete!(instream_exprs, k)
    end

    # This is an implementation of the modelica spec
    # https://specification.modelica.org/maint/3.6/stream-connectors.html
    additional_eqs = Equation[]
    for cset in csets
        n_outer = count(cvert -> cvert.isouter, cset)
        n_inner = length(cset) - n_outer
        if n_inner == 1 && n_outer == 0
            cvert = only(cset)
            stream_var = variable_from_vertex(sys, cvert)::BasicSymbolic
            instream_subs[instream(stream_var)] = stream_var
        elseif n_inner == 2 && n_outer == 0
            cvert1, cvert2 = cset
            stream_var1 = variable_from_vertex(sys, cvert1)::BasicSymbolic
            stream_var2 = variable_from_vertex(sys, cvert2)::BasicSymbolic
            instream_subs[instream(stream_var1)] = stream_var2
            instream_subs[instream(stream_var2)] = stream_var1
        elseif n_inner == 1 && n_outer == 1
            cvert_inner, cvert_outer = cset
            if cvert_inner.isouter
                cvert_inner, cvert_outer = cvert_outer, cvert_inner
            end
            streamvar_inner = variable_from_vertex(sys, cvert_inner)::BasicSymbolic
            streamvar_outer = variable_from_vertex(sys, cvert_outer)::BasicSymbolic
            instream_subs[instream(streamvar_inner)] = instream(streamvar_outer)
            push!(additional_eqs, (streamvar_outer ~ streamvar_inner))
        elseif n_inner == 0 && n_outer == 2
            cvert1, cvert2 = cset
            stream_var1 = variable_from_vertex(sys, cvert1)::BasicSymbolic
            stream_var2 = variable_from_vertex(sys, cvert2)::BasicSymbolic
            push!(additional_eqs, (stream_var1 ~ instream(stream_var2)),
                (stream_var2 ~ instream(stream_var1)))
        else
            # Currently just implements the "else" case for `instream(..)` in the suggested
            # implementation of stream connectors in the Modelica spec v3.6 section 15.2.
            # https://specification.modelica.org/maint/3.6/stream-connectors.html#instream-and-connection-equations
            # We could implement the "if" case using variable bounds? It would be nice to
            # move that metadata to the system (storing it similar to `defaults`).
            outer_cverts = filter(cvert -> cvert.isouter, cset)
            inner_cverts = filter(cvert -> !cvert.isouter, cset)

            outer_streamvars = map(Base.Fix1(variable_from_vertex, sys), outer_cverts)
            inner_streamvars = map(Base.Fix1(variable_from_vertex, sys), inner_cverts)

            outer_flowvars = map(Base.Fix1(get_flowvar, sys), outer_cverts)
            inner_flowvars = map(Base.Fix1(get_flowvar, sys), inner_cverts)

            mask = trues(length(inner_cverts))
            for inner_i in eachindex(inner_cverts)
                # mask out the current variable
                mask[inner_i] = false
                svar = inner_streamvars[inner_i]
                instream_subs[instream(svar)] = term(
                    instream_rt, Val(n_inner - 1), Val(n_outer), inner_flowvars[mask]...,
                    inner_streamvars[mask]..., outer_flowvars..., outer_streamvars...)
                # make sure to reset the mask
                mask[inner_i] = true
            end

            for q in 1:n_outer
                sq = mapreduce(+, inner_flowvars) do fvar
                    max(-fvar, 0)
                end
                sq += mapreduce(+, enumerate(outer_flowvars)) do (outer_i, fvar)
                    outer_i == q && return 0
                    max(fvar, 0)
                end
                # sanity check to make sure it isn't going to codegen a `mapreduce`
                @assert operation(sq) == (+)

                num = mapreduce(+, inner_flowvars, inner_streamvars) do fvar, svar
                    positivemax(-fvar, sq; tol) * svar
                end
                num += mapreduce(
                    +, enumerate(outer_flowvars), outer_streamvars) do (outer_i, fvar), svar
                    outer_i == q && return 0
                    positivemax(fvar, sq; tol) * instream(svar)
                end
                @assert operation(num) == (+)

                den = mapreduce(+, inner_flowvars) do fvar
                    positivemax(-fvar, sq; tol)
                end
                den += mapreduce(+, enumerate(outer_flowvars)) do (outer_i, fvar)
                    outer_i == q && return 0
                    positivemax(fvar, sq; tol)
                end

                push!(additional_eqs, (outer_streamvars[q] ~ num / den))
            end
        end
    end
    return additional_eqs, instream_subs
end

# instream runtime
@generated function _instream_split(::Val{inner_n}, ::Val{outer_n},
        vars::NTuple{N, Any}) where {inner_n, outer_n, N}
    #instream_rt(innerfvs..., innersvs..., outerfvs..., outersvs...)
    ret = Expr(:tuple)
    # mj.c.m_flow
    inner_f = :(Base.@ntuple $inner_n i->vars[i])
    offset = inner_n
    inner_s = :(Base.@ntuple $inner_n i->vars[$offset + i])
    offset += inner_n
    # ck.m_flow
    outer_f = :(Base.@ntuple $outer_n i->vars[$offset + i])
    offset += outer_n
    outer_s = :(Base.@ntuple $outer_n i->vars[$offset + i])
    Expr(:tuple, inner_f, inner_s, outer_f, outer_s)
end

function instream_rt(ins::Val{inner_n}, outs::Val{outer_n},
        vars::Vararg{Any, N}) where {inner_n, outer_n, N}
    @assert N == 2 * (inner_n + outer_n)

    # inner: mj.c.m_flow
    # outer: ck.m_flow
    inner_f, inner_s, outer_f, outer_s = _instream_split(ins, outs, vars)

    T = float(first(inner_f))
    si = zero(T)
    num = den = zero(T)
    for f in inner_f
        si += max(-f, 0)
    end
    for f in outer_f
        si += max(f, 0)
    end
    #for (f, s) in zip(inner_f, inner_s)
    for j in 1:inner_n
        @inbounds f = inner_f[j]
        @inbounds s = inner_s[j]
        num += _positivemax(-f, si) * s
        den += _positivemax(-f, si)
    end
    #for (f, s) in zip(outer_f, outer_s)
    for j in 1:outer_n
        @inbounds f = outer_f[j]
        @inbounds s = outer_s[j]
        num += _positivemax(-f, si) * s
        den += _positivemax(-f, si)
    end
    return num / den
    #=
    si = sum(max(-mj.c.m_flow,0) for j in cat(1,1:i-1, i+1:N)) +
            sum(max(ck.m_flow ,0) for k  in 1:M)

    inStream(mi.c.h_outflow) =
       (sum(positiveMax(-mj.c.m_flow,si)*mj.c.h_outflow)
      +  sum(positiveMax(ck.m_flow,s_i)*inStream(ck.h_outflow)))/
     (sum(positiveMax(-mj.c.m_flow,s_i))
        +  sum(positiveMax(ck.m_flow,s_i)))
                  for j in 1:N and i <> j and mj.c.m_flow.min < 0,
                  for k in 1:M and ck.m_flow.max > 0
    =#
end
SymbolicUtils.promote_symtype(::typeof(instream_rt), ::Vararg) = Real
