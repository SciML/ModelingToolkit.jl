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
    domain_connect(sys1, sys2, syss...)

Adds a domain only connection equation, through and across state equations are not generated.
"""
function domain_connect(sys1, sys2, syss...)
    syss = (sys1, sys2, syss...)
    length(unique(nameof, syss)) == length(syss) || error("connect takes distinct systems!")
    Equation(Connection(:domain), Connection(syss)) # the RHS are connected systems
end

function get_connection_type(s)
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

const CAUSAL_CONNECTION_ERR = """
Only causal variables can be used in a `connect` statement. The first argument must \
be a single output variable and all subsequent variables must be input variables.
"""

function VariableNotOutputError(var)
    ArgumentError("""
    $CAUSAL_CONNECTION_ERR Expected $var to be marked as an output with `[output = true]` \
    in the variable metadata.
    """)
end

function VariableNotInputError(var)
    ArgumentError("""
    $CAUSAL_CONNECTION_ERR Expected $var to be marked an input with `[input = true]` \
    in the variable metadata.
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
    isoutput(var1) || throw(VariableNotOutputError(var1))
    isinput(var2) || throw(VariableNotInputError(var2))
    for var in vars
        isinput(var) || throw(VariableNotInputError(var))
    end
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

function flowvar(sys::AbstractSystem)
    sts = get_unknowns(sys)
    for s in sts
        vtype = get_connection_type(s)
        vtype === Flow && return s
    end
    error("There in no flow variable in $(nameof(sys))")
end

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
            parent_name in outer_connectors
        end
    end
end

struct LazyNamespace
    namespace::Union{Nothing, AbstractSystem}
    sys::Any
end

_getname(::Nothing) = nothing
_getname(sys) = nameof(sys)
Base.copy(l::LazyNamespace) = renamespace(_getname(l.namespace), l.sys)
Base.nameof(l::LazyNamespace) = renamespace(_getname(l.namespace), nameof(l.sys))

struct ConnectionElement
    sys::LazyNamespace
    v::Any
    isouter::Bool
    h::UInt
end
function _hash_impl(sys, v, isouter)
    hashcore = hash(nameof(sys)::Symbol) ⊻ hash(getname(v)::Symbol)
    hashouter = isouter ? hash(true) : hash(false)
    hashcore ⊻ hashouter
end
function ConnectionElement(sys::LazyNamespace, v, isouter::Bool)
    ConnectionElement(sys, v, isouter, _hash_impl(sys, v, isouter))
end
Base.nameof(l::ConnectionElement) = renamespace(nameof(l.sys), getname(l.v))
Base.isequal(l1::ConnectionElement, l2::ConnectionElement) = l1 == l2
function Base.:(==)(l1::ConnectionElement, l2::ConnectionElement)
    l1.isouter == l2.isouter && nameof(l1.sys) == nameof(l2.sys) && isequal(l1.v, l2.v)
end

const _debug_mode = Base.JLOptions().check_bounds == 1

function Base.show(io::IO, c::ConnectionElement)
    @unpack sys, v, isouter = c
    print(io, nameof(sys), ".", v, "::", isouter ? "outer" : "inner")
end

function Base.hash(e::ConnectionElement, salt::UInt)
    if _debug_mode
        @assert e.h === _hash_impl(e.sys, e.v, e.isouter)
    end
    e.h ⊻ salt
end
namespaced_var(l::ConnectionElement) = unknowns(l, l.v)
unknowns(l::ConnectionElement, v) = unknowns(copy(l.sys), v)

function withtrueouter(e::ConnectionElement)
    e.isouter && return e
    # we undo the xor
    newhash = (e.h ⊻ hash(false)) ⊻ hash(true)
    ConnectionElement(e.sys, e.v, true, newhash)
end

struct ConnectionSet
    set::Vector{ConnectionElement} # namespace.sys, var, isouter
end
ConnectionSet() = ConnectionSet(ConnectionElement[])
Base.copy(c::ConnectionSet) = ConnectionSet(copy(c.set))
Base.:(==)(a::ConnectionSet, b::ConnectionSet) = a.set == b.set
Base.sort(a::ConnectionSet) = ConnectionSet(sort(a.set, by = string))

function Base.show(io::IO, c::ConnectionSet)
    print(io, "<")
    for i in 1:(length(c.set) - 1)
        @unpack sys, v, isouter = c.set[i]
        print(io, nameof(sys), ".", v, "::", isouter ? "outer" : "inner", ", ")
    end
    @unpack sys, v, isouter = last(c.set)
    print(io, nameof(sys), ".", v, "::", isouter ? "outer" : "inner", ">")
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
    val = getmetadata(sys, FrameOrientation, nothing)
    if val === nothing
        error("System $(sys.name) does not have an orientation object.")
    end
end

"""
    $(TYPEDSIGNATURES)

Populate `connectionsets` with connections between the connectors `ss`, all of which are
namespaced by `namespace`.

# Keyword Arguments
- `ignored_connects`: A tuple of the systems and variables for which connections should be
  ignored. Of the format returned from `as_hierarchy`.
- `namespaced_ignored_systems`: The `from_hierarchy` versions of entries in
  `ignored_connects[1]`, purely to avoid unnecessary recomputation.
"""
function connection2set!(connectionsets, namespace, ss, isouter;
        ignored_systems = HierarchySystemT[], ignored_variables = HierarchyVariableT[])
    ns_ignored_systems = from_hierarchy.(ignored_systems)
    ns_ignored_variables = from_hierarchy.(ignored_variables)
    # ignore specified systems
    ss = filter(ss) do s
        !any(x -> nameof(x) == nameof(s), ns_ignored_systems)
    end
    # `ignored_variables` for each `s` in `ss`
    corresponding_ignored_variables = map(
        Base.Fix2(ignored_systems_for_subsystem, ignored_variables), ss)
    corresponding_namespaced_ignored_variables = map(
        Broadcast.BroadcastFunction(from_hierarchy), corresponding_ignored_variables)

    regular_ss = []
    domain_ss = nothing
    for s in ss
        if is_domain_connector(s)
            if domain_ss === nothing
                domain_ss = s
            else
                names = join(map(string ∘ nameof, ss), ",")
                error("connect($names) contains multiple source domain connectors. There can only be one!")
            end
        else
            push!(regular_ss, s)
        end
    end
    T = ConnectionElement
    @assert !isempty(regular_ss)
    ss = regular_ss
    # domain connections don't generate any equations
    if domain_ss !== nothing
        cset = ConnectionElement[]
        dv = only(unknowns(domain_ss))
        for (i, s) in enumerate(ss)
            sts = unknowns(s)
            io = isouter(s)
            _ignored_variables = corresponding_ignored_variables[i]
            _namespaced_ignored_variables = corresponding_namespaced_ignored_variables[i]
            for v in sts
                vtype = get_connection_type(v)
                (vtype === Flow && isequal(v, dv)) || continue
                any(isequal(v), _namespaced_ignored_variables) && continue
                push!(cset, T(LazyNamespace(namespace, domain_ss), dv, false))
                push!(cset, T(LazyNamespace(namespace, s), v, io))
            end
        end
        @assert length(cset) > 0
        push!(connectionsets, ConnectionSet(cset))
        return connectionsets
    end
    s1 = first(ss)
    sts1v = unknowns(s1)
    if isframe(s1) # Multibody
        O = ori(s1)
        orientation_vars = Symbolics.unwrap.(collect(vec(O.R)))
        sts1v = [sts1v; orientation_vars]
    end
    sts1 = Set(sts1v)
    num_unknowns = length(sts1)

    # we don't filter here because `csets` should include the full set of unknowns.
    # not all of `ss` will have the same (or any) variables filtered so the ones
    # that aren't should still go in the right cset. Since `sts1` is only used for
    # validating that all systems being connected are of the same type, it has
    # unfiltered entries.
    csets = [T[] for _ in 1:num_unknowns] # Add 9 orientation variables if connection is between multibody frames
    for (i, s) in enumerate(ss)
        unknown_vars = unknowns(s)
        if isframe(s) # Multibody
            O = ori(s)
            orientation_vars = Symbolics.unwrap.(vec(O.R))
            unknown_vars = [unknown_vars; orientation_vars]
        end
        i != 1 && ((num_unknowns == length(unknown_vars) &&
          all(Base.Fix2(in, sts1), unknown_vars)) ||
         connection_error(ss))
        io = isouter(s)
        # don't `filter!` here so that `j` points to the correct cset regardless of
        # which variables are filtered.
        for (j, v) in enumerate(unknown_vars)
            any(isequal(v), corresponding_namespaced_ignored_variables[i]) && continue
            push!(csets[j], T(LazyNamespace(namespace, s), v, io))
        end
    end
    for cset in csets
        v = first(cset).v
        vtype = get_connection_type(v)
        if domain_ss !== nothing && vtype === Flow &&
           (dv = only(unknowns(domain_ss)); isequal(v, dv))
            push!(cset, T(LazyNamespace(namespace, domain_ss), dv, false))
        end
        for k in 2:length(cset)
            vtype === get_connection_type(cset[k].v) || connection_error(ss)
        end
        push!(connectionsets, ConnectionSet(cset))
    end
end

function generate_connection_set(
        sys::AbstractSystem; scalarize = false)
    connectionsets = ConnectionSet[]
    domain_csets = ConnectionSet[]
    sys = generate_connection_set!(
        connectionsets, domain_csets, sys, scalarize, nothing, ignored_connections(sys))
    csets = merge(connectionsets)
    domain_csets = merge([csets; domain_csets], true)

    sys, (csets, domain_csets)
end

"""
    $(TYPEDSIGNATURES)

For a list of `systems` in a connect equation, return the subset of it to ignore (as a
list of hierarchical systems) based on `ignored_system_aps`, the analysis points to be
ignored. All analysis points in `ignored_system_aps` must contain systems (connectors)
as their input/outputs.
"""
function systems_to_ignore(ignored_system_aps::Vector{HierarchyAnalysisPointT},
        systems::Union{Vector{S}, Tuple{Vararg{S}}}) where {S <: AbstractSystem}
    to_ignore = HierarchySystemT[]
    for ap in ignored_system_aps
        # if `systems` contains the input of the AP, ignore any outputs of the AP present in it.
        isys_hierarchy = HierarchySystemT([ap[1].input; @view ap[2:end]])
        isys = from_hierarchy(isys_hierarchy)
        any(x -> nameof(x) == nameof(isys), systems) || continue

        for outsys in ap[1].outputs
            osys_hierarchy = HierarchySystemT([outsys; @view ap[2:end]])
            osys = from_hierarchy(osys_hierarchy)
            any(x -> nameof(x) == nameof(osys), systems) || continue
            push!(to_ignore, HierarchySystemT(osys_hierarchy))
        end
    end

    return to_ignore
end

"""
    $(TYPEDSIGNATURES)

For a list of `systems` in a connect equation, return the subset of their variables to
ignore (as a list of hierarchical variables) based on `ignored_system_aps`, the analysis
points to be ignored. All analysis points in `ignored_system_aps` must contain variables
as their input/outputs.
"""
function variables_to_ignore(ignored_variable_aps::Vector{HierarchyAnalysisPointT},
        systems::Union{Vector{S}, Tuple{Vararg{S}}}) where {S <: AbstractSystem}
    to_ignore = HierarchyVariableT[]
    for ap in ignored_variable_aps
        ivar_hierarchy = HierarchyVariableT([ap[1].input; @view ap[2:end]])
        ivar = from_hierarchy(ivar_hierarchy)
        any(x -> any(isequal(ivar), renamespace.((x,), unknowns(x))), systems) || continue

        for outvar in ap[1].outputs
            ovar_hierarchy = HierarchyVariableT([as_hierarchy(outvar); @view ap[2:end]])
            ovar = from_hierarchy(ovar_hierarchy)
            any(x -> any(isequal(ovar), renamespace.((x,), unknowns(x))), systems) ||
                continue
            push!(to_ignore, HierarchyVariableT(ovar_hierarchy))
        end
    end
    return to_ignore
end

"""
    $(TYPEDSIGNATURES)

For a list of variables `vars` in a connect equation, return the subset of them ignore
(as a list of symbolic variables) based on `ignored_system_aps`, the analysis points to
be ignored. All analysis points in `ignored_system_aps` must contain variables as their
input/outputs.
"""
function variables_to_ignore(ignored_variable_aps::Vector{HierarchyAnalysisPointT},
        vars::Union{Vector{S}, Tuple{Vararg{S}}}) where {S <: BasicSymbolic}
    to_ignore = eltype(vars)[]
    for ap in ignored_variable_aps
        ivar_hierarchy = HierarchyVariableT([ap[1].input; @view ap[2:end]])
        ivar = from_hierarchy(ivar_hierarchy)
        any(isequal(ivar), vars) || continue

        for outvar in ap[1].outputs
            ovar_hierarchy = HierarchyVariableT([outvar; @view ap[2:end]])
            ovar = from_hierarchy(ovar_hierarchy)
            any(isequal(ovar), vars) || continue
            push!(to_ignore, ovar)
        end
    end

    return to_ignore
end

"""
    $(TYPEDSIGNATURES)

Generate connection sets from `connect` equations.

# Arguments

- `connectionsets` is the list of connection sets to be populated by recursively
  descending `sys`.
- `domain_csets` is the list of connection sets for domain connections.
- `sys` is the system whose equations are to be searched.
- `namespace` is a system representing the namespace in which `sys` exists, or `nothing`
  for no namespace (if `sys` is top-level).
"""
function generate_connection_set!(connectionsets, domain_csets,
        sys::AbstractSystem, scalarize, namespace = nothing,
        ignored_connects = (HierarchyAnalysisPointT[], HierarchyAnalysisPointT[]))
    subsys = get_systems(sys)
    ignored_system_aps, ignored_variable_aps = ignored_connects

    isouter = generate_isouter(sys)
    eqs′ = get_eqs(sys)
    eqs = Equation[]

    cts = [] # connections
    domain_cts = [] # connections
    extra_unknowns = []
    for eq in eqs′
        lhs = eq.lhs
        rhs = eq.rhs

        # causal variable connections will be expanded before we get here,
        # but this guard is useful for `n_expanded_connection_equations`.
        is_causal_variable_connection(rhs) && continue
        if lhs isa Connection && get_systems(lhs) === :domain
            connected_systems = get_systems(rhs)
            connection2set!(domain_csets, namespace, connected_systems, isouter;
                ignored_systems = systems_to_ignore(
                    ignored_system_aps, connected_systems),
                ignored_variables = variables_to_ignore(
                    ignored_variable_aps, connected_systems))
        elseif isconnection(rhs)
            push!(cts, get_systems(rhs))
        else
            # split connections and equations
            if eq.lhs isa AbstractArray || eq.rhs isa AbstractArray
                append!(eqs, Symbolics.scalarize(eq))
            else
                push!(eqs, eq)
            end
        end
    end

    # all connectors are eventually inside connectors.
    T = ConnectionElement
    # only generate connection sets for systems that are not ignored
    for s in subsys
        isconnector(s) || continue
        is_domain_connector(s) && continue
        for v in unknowns(s)
            Flow === get_connection_type(v) || continue
            push!(connectionsets, ConnectionSet([T(LazyNamespace(namespace, s), v, false)]))
        end
    end

    for ct in cts
        connection2set!(connectionsets, namespace, ct, isouter;
            ignored_systems = systems_to_ignore(ignored_system_aps, ct),
            ignored_variables = variables_to_ignore(ignored_variable_aps, ct))
    end

    # pre order traversal
    if !isempty(extra_unknowns)
        @set! sys.unknowns = [get_unknowns(sys); extra_unknowns]
    end
    @set! sys.systems = map(
        s -> generate_connection_set!(connectionsets, domain_csets, s,
            scalarize, renamespace(namespace, s),
            ignored_systems_for_subsystem.((s,), ignored_connects)),
        subsys)
    @set! sys.eqs = eqs
end

"""
    $(TYPEDSIGNATURES)

Given a subsystem `subsys` of a parent system and a list of systems (variables) to be
ignored by `generate_connection_set!` (`expand_variable_connections`), filter
`ignored_systems` to only include those present in the subtree of `subsys` and update
their hierarchy to not include `subsys`.
"""
function ignored_systems_for_subsystem(
        subsys::AbstractSystem, ignored_systems::Vector{<:Union{
            HierarchyT, HierarchyAnalysisPointT}})
    result = eltype(ignored_systems)[]
    # in case `subsys` is namespaced, get its hierarchy and compare suffixes
    # instead of the just the last element
    suffix = reverse!(namespace_hierarchy(nameof(subsys)))
    N = length(suffix)
    for igsys in ignored_systems
        if length(igsys) > N && igsys[(end - N + 1):end] == suffix
            push!(result, copy(igsys))
            for i in 1:N
                pop!(result[end])
            end
        end
    end
    return result
end

function Base.merge(csets::AbstractVector{<:ConnectionSet}, allouter = false)
    ele2idx = Dict{ConnectionElement, Int}()
    idx2ele = ConnectionElement[]
    union_find = IntDisjointSets(0)
    prev_id = Ref(-1)
    for cset in csets, (j, s) in enumerate(cset.set)
        v = allouter ? withtrueouter(s) : s
        id = let ele2idx = ele2idx, idx2ele = idx2ele
            get!(ele2idx, v) do
                push!(idx2ele, v)
                id = length(idx2ele)
                id′ = push!(union_find)
                @assert id == id′
                id
            end
        end
        # isequal might not be equal? lol
        if v.sys.namespace !== nothing
            idx2ele[id] = v
        end
        if j > 1
            union!(union_find, prev_id[], id)
        end
        prev_id[] = id
    end
    id2set = Dict{Int, Int}()
    merged_set = ConnectionSet[]
    for (id, ele) in enumerate(idx2ele)
        rid = find_root!(union_find, id)
        set_idx = get!(id2set, rid) do
            set = ConnectionSet()
            push!(merged_set, set)
            length(merged_set)
        end
        push!(merged_set[set_idx].set, ele)
    end
    merged_set
end

function generate_connection_equations_and_stream_connections(csets::AbstractVector{
        <:ConnectionSet,
})
    eqs = Equation[]
    stream_connections = ConnectionSet[]

    for cset in csets
        v = cset.set[1].v
        v = getparent(v, v)
        vtype = get_connection_type(v)
        if vtype === Stream
            push!(stream_connections, cset)
        elseif vtype === Flow
            rhs = 0
            for ele in cset.set
                v = namespaced_var(ele)
                rhs += ele.isouter ? -v : v
            end
            push!(eqs, 0 ~ rhs)
        else # Equality
            base = namespaced_var(cset.set[1])
            for i in 2:length(cset.set)
                v = namespaced_var(cset.set[i])
                push!(eqs, base ~ v)
            end
        end
    end
    eqs, stream_connections
end

function domain_defaults(sys, domain_csets)
    def = Dict()
    for c in domain_csets
        cset = c.set
        idx = findfirst(s -> is_domain_connector(s.sys.sys), cset)
        idx === nothing && continue
        s = cset[idx]
        root = s.sys
        s_def = defaults(root.sys)
        for (j, m) in enumerate(cset)
            if j == idx
                continue
            elseif is_domain_connector(m.sys.sys)
                error("Domain sources $(nameof(root)) and $(nameof(m)) are connected!")
            else
                ns_s_def = Dict(unknowns(m.sys.sys, n) => n for (n, v) in s_def)
                for p in parameters(m.sys.namespace)
                    d_p = get(ns_s_def, p, nothing)
                    if d_p !== nothing
                        def[parameters(m.sys.namespace, p)] = parameters(s.sys.namespace,
                            parameters(s.sys.sys,
                                d_p))
                    end
                end
            end
        end
    end
    def
end

"""
    $(TYPEDSIGNATURES)

Recursively descend through the hierarchy of `sys` and expand all connection equations
of causal variables. Return the modified system.
"""
function expand_variable_connections(sys::AbstractSystem; ignored_variables = nothing)
    if ignored_variables === nothing
        ignored_variables = ignored_connections(sys)[2]
    end
    eqs = copy(get_eqs(sys))
    valid_idxs = trues(length(eqs))
    additional_eqs = Equation[]

    for (i, eq) in enumerate(eqs)
        eq.lhs isa Connection || continue
        connection = eq.rhs
        elements = get_systems(connection)
        is_causal_variable_connection(connection) || continue

        valid_idxs[i] = false
        elements = map(x -> x.var, elements)
        to_ignore = variables_to_ignore(ignored_variables, elements)
        elements = setdiff(elements, to_ignore)
        outvar = first(elements)
        for invar in Iterators.drop(elements, 1)
            push!(additional_eqs, outvar ~ invar)
        end
    end
    eqs = [eqs[valid_idxs]; additional_eqs]
    subsystems = map(get_systems(sys)) do subsys
        expand_variable_connections(subsys;
            ignored_variables = ignored_systems_for_subsystem(subsys, ignored_variables))
    end
    @set! sys.eqs = eqs
    @set! sys.systems = subsystems
    return sys
end

"""
    function expand_connections(sys::AbstractSystem)

Given a hierarchical system with [`connect`](@ref) equations, expand the connection
equations and return the new system.
"""
function expand_connections(sys::AbstractSystem;
        debug = false, tol = 1e-10, scalarize = true)
    sys = remove_analysis_points(sys)
    sys = expand_variable_connections(sys)
    sys, (csets, domain_csets) = generate_connection_set(sys; scalarize)
    ceqs, instream_csets = generate_connection_equations_and_stream_connections(csets)
    _sys = expand_instream(instream_csets, sys; debug = debug, tol = tol)
    sys = flatten(sys, true)
    @set! sys.eqs = [equations(_sys); ceqs]
    d_defs = domain_defaults(sys, domain_csets)
    @set! sys.defaults = merge(get_defaults(sys), d_defs)
end

function unnamespace(root, namespace)
    root === nothing && return namespace
    root = string(root)
    namespace = string(namespace)
    if length(namespace) > length(root)
        @assert root == namespace[1:length(root)]
        Symbol(namespace[nextind(namespace, length(root)):end])
    else
        @assert root == namespace
        nothing
    end
end

function expand_instream(csets::AbstractVector{<:ConnectionSet}, sys::AbstractSystem,
        namespace = nothing, prevnamespace = nothing; debug = false,
        tol = 1e-8)
    subsys = get_systems(sys)
    # post order traversal
    @set! sys.systems = map(
        s -> expand_instream(csets, s,
            renamespace(namespace, nameof(s)),
            namespace; debug, tol),
        subsys)
    subsys = get_systems(sys)

    if debug
        @info "Expanding" namespace
    end

    sub = Dict()
    eqs = Equation[]
    instream_eqs = Equation[]
    instream_exprs = Set()
    for s in subsys
        for eq in get_eqs(s)
            eq = namespace_equation(eq, s)
            if collect_instream!(instream_exprs, eq)
                push!(instream_eqs, eq)
            else
                push!(eqs, eq)
            end
        end
    end

    if !isempty(instream_exprs)
        # map from a namespaced stream variable to a ConnectionSet
        expr_cset = Dict()
        for cset in csets
            crep = first(cset.set)
            current = namespace == _getname(crep.sys.namespace)
            for v in cset.set
                if (current || !v.isouter)
                    expr_cset[namespaced_var(v)] = cset.set
                end
            end
        end
    end

    for ex in instream_exprs
        ns_sv = only(arguments(ex))
        full_name_sv = renamespace(namespace, ns_sv)
        cset = get(expr_cset, full_name_sv, nothing)
        cset === nothing && error("$ns_sv is not a variable inside stream connectors")
        idx_in_set, sv = get_cset_sv(full_name_sv, cset)

        n_inners = n_outers = 0
        for (i, e) in enumerate(cset)
            if e.isouter
                n_outers += 1
            else
                n_inners += 1
            end
        end
        if debug
            @info "Expanding at [$idx_in_set]" ex ConnectionSet(cset) n_inners n_outers
        end
        if n_inners == 1 && n_outers == 0
            sub[ex] = sv
        elseif n_inners == 2 && n_outers == 0
            other = idx_in_set == 1 ? 2 : 1
            sub[ex] = get_current_var(namespace, cset[other], sv)
        elseif n_inners == 1 && n_outers == 1
            if !cset[idx_in_set].isouter
                other = idx_in_set == 1 ? 2 : 1
                outerstream = get_current_var(namespace, cset[other], sv)
                sub[ex] = instream(outerstream)
            end
        else
            if !cset[idx_in_set].isouter
                fv = flowvar(first(cset).sys.sys)
                # mj.c.m_flow
                innerfvs = [get_current_var(namespace, s, fv)
                            for (j, s) in enumerate(cset) if j != idx_in_set && !s.isouter]
                innersvs = [get_current_var(namespace, s, sv)
                            for (j, s) in enumerate(cset) if j != idx_in_set && !s.isouter]
                # ck.m_flow
                outerfvs = [get_current_var(namespace, s, fv) for s in cset if s.isouter]
                outersvs = [get_current_var(namespace, s, sv) for s in cset if s.isouter]

                sub[ex] = term(instream_rt, Val(length(innerfvs)), Val(length(outerfvs)),
                    innerfvs..., innersvs..., outerfvs..., outersvs...)
            end
        end
    end

    # additional equations
    additional_eqs = Equation[]
    csets = filter(cset -> any(e -> _getname(e.sys.namespace) === namespace, cset.set),
        csets)
    for cset′ in csets
        cset = cset′.set
        connectors = Vector{Any}(undef, length(cset))
        n_inners = n_outers = 0
        for (i, e) in enumerate(cset)
            connectors[i] = e.sys.sys
            if e.isouter
                n_outers += 1
            else
                n_inners += 1
            end
        end
        iszero(n_outers) && continue
        connector_representative = first(cset).sys.sys
        fv = flowvar(connector_representative)
        sv = first(cset).v
        vtype = get_connection_type(sv)
        vtype === Stream || continue
        if n_inners == 1 && n_outers == 1
            push!(additional_eqs,
                unknowns(cset[1].sys.sys, sv) ~ unknowns(cset[2].sys.sys, sv))
        elseif n_inners == 0 && n_outers == 2
            # we don't expand `instream` in this case.
            v1 = unknowns(cset[1].sys.sys, sv)
            v2 = unknowns(cset[2].sys.sys, sv)
            push!(additional_eqs, v1 ~ instream(v2))
            push!(additional_eqs, v2 ~ instream(v1))
        else
            sq = 0
            s_inners = (s for s in cset if !s.isouter)
            s_outers = (s for s in cset if s.isouter)
            for (q, oscq) in enumerate(s_outers)
                sq += sum(s -> max(-unknowns(s, fv), 0), s_inners, init = 0)
                for (k, s) in enumerate(s_outers)
                    k == q && continue
                    f = unknowns(s.sys.sys, fv)
                    sq += max(f, 0)
                end

                num = 0
                den = 0
                for s in s_inners
                    f = unknowns(s.sys.sys, fv)
                    tmp = positivemax(-f, sq; tol = tol)
                    den += tmp
                    num += tmp * unknowns(s.sys.sys, sv)
                end
                for (k, s) in enumerate(s_outers)
                    k == q && continue
                    f = unknowns(s.sys.sys, fv)
                    tmp = positivemax(f, sq; tol = tol)
                    den += tmp
                    num += tmp * instream(unknowns(s.sys.sys, sv))
                end
                push!(additional_eqs, unknowns(oscq.sys.sys, sv) ~ num / den)
            end
        end
    end

    subed_eqs = substitute(instream_eqs, sub)
    if debug && !(isempty(csets) && isempty(additional_eqs) && isempty(instream_eqs))
        println("======================================")
        @info "Additional equations" csets
        display(additional_eqs)
        println("======================================")
        println("Substitutions")
        display(sub)
        println("======================================")
        println("Substituted equations")
        foreach(i -> println(instream_eqs[i] => subed_eqs[i]), eachindex(subed_eqs))
        println("======================================")
    end

    @set! sys.systems = []
    @set! sys.eqs = [get_eqs(sys); eqs; subed_eqs; additional_eqs]
    sys
end

function get_current_var(namespace, cele, sv)
    unknowns(
        renamespace(unnamespace(namespace, _getname(cele.sys.namespace)),
            cele.sys.sys),
        sv)
end

function get_cset_sv(full_name_sv, cset)
    for (idx_in_set, v) in enumerate(cset)
        if isequal(namespaced_var(v), full_name_sv)
            return idx_in_set, v.v
        end
    end
    error("$ns_sv is not a variable inside stream connectors")
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
