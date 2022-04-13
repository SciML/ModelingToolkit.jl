get_connection_type(s) = getmetadata(unwrap(s), VariableConnectType, Equality)

function with_connector_type(expr)
    @assert expr isa Expr && (expr.head == :function || (expr.head == :(=) &&
                                       expr.args[1] isa Expr &&
                                       expr.args[1].head == :call))

    sig = expr.args[1]
    body = expr.args[2]

    fname = sig.args[1]
    args = sig.args[2:end]

    quote
        function $fname($(args...))
            function f()
                $body
            end
            res = f()
            $isdefined(res, :connector_type) && $getfield(res, :connector_type) === nothing ? $Setfield.@set!(res.connector_type = $connector_type(res)) : res
        end
    end
end

macro connector(expr)
    esc(with_connector_type(expr))
end

abstract type AbstractConnectorType end
struct StreamConnector <: AbstractConnectorType end
struct RegularConnector <: AbstractConnectorType end

function connector_type(sys::AbstractSystem)
    sts = get_states(sys)
    #TODO: check the criteria for stream connectors
    n_stream = 0
    n_flow = 0
    for s in sts
        vtype = get_connection_type(s)
        if vtype === Stream
            isarray(s) && error("Array stream variables are not supported. Got $s.")
            n_stream += 1
        end
        vtype === Flow && (n_flow += 1)
    end
    (n_stream > 0 && n_flow > 1) && error("There are multiple flow variables in $(nameof(sys))!")
    n_stream > 0 ? StreamConnector() : RegularConnector()
end

Base.@kwdef struct Connection
    inners = nothing
    outers = nothing
end

# everything is inner by default until we expand the connections
Connection(syss) = Connection(inners=syss)
get_systems(c::Connection) = c.inners
function Base.in(e::Symbol, c::Connection)
    (c.inners !== nothing && any(k->nameof(k) === e, c.inners)) ||
    (c.outers !== nothing && any(k->nameof(k) === e, c.outers))
end

function renamespace(sym::Symbol, connection::Connection)
    inners = connection.inners === nothing ? [] : renamespace.(sym, connection.inners)
    if connection.outers !== nothing
        for o in connection.outers
            push!(inners, renamespace(sym, o))
        end
    end
    Connection(;inners=inners)
end

const EMPTY_VEC = []

function Base.show(io::IO, ::MIME"text/plain", c::Connection)
    # It is a bit unfortunate that the display of an array of `Equation`s won't
    # call this.
    @unpack outers, inners = c
    if outers === nothing && inners === nothing
        print(io, "<Connection>")
    else
        syss = Iterators.flatten((something(inners, EMPTY_VEC), something(outers, EMPTY_VEC)))
        splitting_idx = length(inners)
        sys_str = join((string(nameof(s)) * (i <= splitting_idx ? ("::inner") : ("::outer")) for (i, s) in enumerate(syss)), ", ")
        print(io, "<", sys_str, ">")
    end
end

# symbolic `connect`
function connect(sys1::AbstractSystem, sys2::AbstractSystem, syss::AbstractSystem...)
    syss = (sys1, sys2, syss...)
    length(unique(nameof, syss)) == length(syss) || error("connect takes distinct systems!")
    Equation(Connection(), Connection(syss)) # the RHS are connected systems
end

# the actual `connect`.
function connect(c::Connection; check=true)
    @unpack inners, outers = c

    cnts = Iterators.flatten((inners, outers))
    fs, ss = Iterators.peel(cnts)
    splitting_idx = length(inners) # anything after the splitting_idx is outer.
    first_sts = get_states(fs)
    first_sts_set = Set(getname.(first_sts))
    for sys in ss
        current_sts = getname.(get_states(sys))
        Set(current_sts) == first_sts_set || error("$(nameof(sys)) ($current_sts) doesn't match the connection type of $(nameof(fs)) ($first_sts).")
    end

    seen = Set()
    ceqs = Equation[]
    for s in first_sts
        name = getname(s)
        fix_val = getproperty(fs, name) # representative
        fix_val in seen && continue
        push!(seen, fix_val)

        vtype = get_connection_type(fix_val)
        vtype === Stream && continue

        if vtype === Flow
            for j in eachindex(fix_val)
                rhs = 0
                for (i, c) in enumerate(cnts)
                    isinner = i <= splitting_idx
                    var = getproperty(c, name)
                    rhs += isinner ? var[j] : -var[j]
                end
                push!(ceqs, 0 ~ rhs)
            end
        else # Equality
            for c in ss
                var = getproperty(c, name)
                for (i, v) in enumerate(var)
                    push!(ceqs, fix_val[i] ~ v)
                end
            end
        end
    end

    return ceqs
end

instream(a) = term(instream, unwrap(a), type=symtype(a))
SymbolicUtils.promote_symtype(::typeof(instream), _) = Real

isconnector(s::AbstractSystem) = has_connector_type(s) && get_connector_type(s) !== nothing
isstreamconnector(s::AbstractSystem) = isconnector(s) && get_connector_type(s) isa StreamConnector
isstreamconnection(c::Connection) = any(isstreamconnector, c.inners) || any(isstreamconnector, c.outers)

function print_with_indent(n, x)
    print(" " ^ n)
    show(stdout, MIME"text/plain"(), x)
    println()
end

function split_sys_var(var)
    var_name = string(getname(var))
    sidx = findlast(isequal('₊'), var_name)
    sidx === nothing && error("$var is not a namespaced variable")
    connector_name = Symbol(var_name[1:prevind(var_name, sidx)])
    streamvar_name = Symbol(var_name[nextind(var_name, sidx):end])
    connector_name, streamvar_name
end

function flowvar(sys::AbstractSystem)
    sts = get_states(sys)
    for s in sts
        vtype = get_connection_type(s)
        vtype === Flow && return s
    end
    error("There in no flow variable in $(nameof(sys))")
end

collect_instream!(set, eq::Equation) = collect_instream!(set, eq.lhs) | collect_instream!(set, eq.rhs)

function collect_instream!(set, expr, occurs=false)
    istree(expr) || return occurs
    op = operation(expr)
    op === instream && (push!(set, expr); occurs = true)
    for a in SymbolicUtils.unsorted_arguments(expr)
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
            (si/eps)^2*(3-2* si/eps)
        else
            zero(T)
        end
    end
    alpha * max(m, 0) + (1-alpha)*eps
end
@register _positivemax(m, tol)
positivemax(m, ::Any; tol=nothing) = _positivemax(m, tol)
mydiv(num, den) = if den == 0
    error()
else
    num / den
end
@register mydiv(n, d)

function generate_isouter(sys::AbstractSystem)
    outer_connectors = Symbol[]
    for s in get_systems(sys)
        n = nameof(s)
        isconnector(s) && push!(outer_connectors, n)
    end
    let outer_connectors=outer_connectors
        function isouter(sys)::Bool
            s = string(nameof(sys))
            isconnector(sys) || error("$s is not a connector!")
            idx = findfirst(isequal('₊'), s)
            parent_name = Symbol(idx === nothing ? s : s[1:prevind(s, idx)])
            parent_name in outer_connectors
        end
    end
end

struct LazyNamespace
    namespace::Union{Nothing,Symbol}
    sys
end

Base.copy(l::LazyNamespace) = renamespace(l.namespace, l.sys)
Base.nameof(l::LazyNamespace) = renamespace(l.namespace, nameof(l.sys))

struct ConnectionElement
    sys::LazyNamespace
    v
    isouter::Bool
end
Base.hash(l::ConnectionElement, salt::UInt) = hash(nameof(l.sys)) ⊻ hash(l.v) ⊻ hash(l.isouter) ⊻ salt
Base.isequal(l1::ConnectionElement, l2::ConnectionElement) = l1 == l2
Base.:(==)(l1::ConnectionElement, l2::ConnectionElement) = nameof(l1.sys) == nameof(l2.sys) && isequal(l1.v, l2.v) && l1.isouter == l2.isouter
namespaced_var(l::ConnectionElement) = states(l, l.v)
states(l::ConnectionElement, v) = states(copy(l.sys), v)

struct ConnectionSet
    set::Vector{ConnectionElement} # namespace.sys, var, isouter
end

function Base.show(io::IO, c::ConnectionSet)
    print(io, "<")
    for i in 1:length(c.set)-1
        @unpack sys, v, isouter = c.set[i]
        print(io, nameof(sys), ".", v, "::", isouter ? "outer" : "inner", ", ")
    end
    @unpack sys, v, isouter = last(c.set)
    print(io, nameof(sys), ".", v, "::", isouter ? "outer" : "inner", ">")
end

@noinline connection_error(ss) = error("Different types of connectors are in one conenction statement: <$(map(nameof, ss))>")

function connection2set!(connectionsets, namespace, ss, isouter)
    nn = map(nameof, ss)
    sts1 = Set(states(first(ss)))
    T = ConnectionElement
    csets = [T[] for _ in 1:length(sts1)]
    for (i, s) in enumerate(ss)
        sts = states(s)
        i != 1 && ((length(sts1) == length(sts) && all(Base.Fix2(in, sts1), sts)) || connection_error(ss))
        io = isouter(s)
        for (j, v) in enumerate(sts)
            push!(csets[j], T(LazyNamespace(namespace, s), v, io))
        end
    end
    for cset in csets
        vtype = get_connection_type(first(cset).v)
        for k in 2:length(cset)
            vtype === get_connection_type(cset[k].v) || connection_error(ss)
        end
        push!(connectionsets, ConnectionSet(cset))
    end
end

function generate_connection_set(sys::AbstractSystem)
    connectionsets = ConnectionSet[]
    sys = generate_connection_set!(connectionsets, sys)
    sys, merge(connectionsets)
end

function generate_connection_set!(connectionsets, sys::AbstractSystem, namespace=nothing)
    @show namespace

    subsys = get_systems(sys)

    isouter = generate_isouter(sys)
    eqs′ = get_eqs(sys)
    eqs = Equation[]

    cts = [] # connections
    for eq in eqs′
        if eq.lhs isa Connection
            push!(cts, get_systems(eq.rhs))
        else
            push!(eqs, eq) # split connections and equations
        end
    end

    if namespace !== nothing
        # Except for the top level, all connectors are eventually inside
        # connectors.
        T = ConnectionElement
        for s in subsys
            isconnector(s) || continue
            for v in states(s)
                Flow === get_connection_type(v) || continue
                push!(connectionsets, ConnectionSet([T(LazyNamespace(namespace, s), v, false)]))
            end
        end
    end

    for ct in cts
        connection2set!(connectionsets, namespace, ct, isouter)
    end

    # pre order traversal
    @set! sys.systems = map(s->generate_connection_set!(connectionsets, s, renamespace(namespace, nameof(s))), subsys)
    @set! sys.eqs = eqs
end

function Base.merge(csets::AbstractVector{<:ConnectionSet})
    mcsets = ConnectionSet[]
    # FIXME: this is O(m n^3)
    for cset in csets
        idx = findfirst(mcset->any(s->any(z->z == s, cset.set), mcset.set), mcsets)
        if idx === nothing
            push!(mcsets, cset)
        else
            union!(mcsets[idx].set, cset.set)
        end
    end
    mcsets
end

function generate_connection_equations_and_stream_connections(csets::AbstractVector{<:ConnectionSet})
    eqs = Equation[]
    stream_connections = ConnectionSet[]

    for cset in csets
        v = cset.set[1].v
        if hasmetadata(v, Symbolics.GetindexParent)
            v = getparent(v)
        end
        vtype = get_connection_type(v)
        if vtype === Stream
            push!(stream_connections, cset)
            continue
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

function expand_connections(sys::AbstractSystem; debug=false, tol=1e-10)
    sys, csets = generate_connection_set(sys)
    ceqs, instream_csets = generate_connection_equations_and_stream_connections(csets)
    additional_eqs = Equation[]
    _sys = expand_instream2(instream_csets, sys; debug=debug, tol=tol)
    Main._a[] = _sys
    sys = flatten(sys)
    @set! sys.eqs = [equations(_sys); ceqs; additional_eqs]
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

function expand_instream2(csets::AbstractVector{<:ConnectionSet}, sys::AbstractSystem, namespace=nothing, prevnamespace=nothing; debug=false, tol=1e-8)
    subsys = get_systems(sys)
    # no connectors if there are no subsystems
    #isempty(subsys) && return sys
    # post order traversal
    @set! sys.systems = map(s->expand_instream2(csets, s, renamespace(namespace, nameof(s)), namespace; debug, tol), subsys)
    subsys = get_systems(sys)
    @info "Expanding" namespace

    sub = Dict()
    eqs = Equation[]
    instream_eqs = Equation[]
    instream_exprs = Set()
    for s in subsys
        seqs = map(Base.Fix2(namespace_equation, s), get_eqs(s))
        for eq in seqs
            if collect_instream!(instream_exprs, eq)
                push!(instream_eqs, eq)
            else
                push!(eqs, eq)
            end
        end

    end

    for ex in instream_exprs
        cset, idx_in_set, sv = get_cset_sv(namespace, ex, csets)

        n_inners = n_outers = 0
        for (i, e) in enumerate(cset)
            if e.isouter
                n_outers += 1
            else
                n_inners += 1
            end
        end
        @show ex idx_in_set ConnectionSet(cset)
        @show n_inners, n_outers
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
                innerfvs = [get_current_var(namespace, s, fv) for (j, s) in enumerate(cset) if j != idx_in_set && !s.isouter]
                innersvs = [get_current_var(namespace, s, sv) for (j, s) in enumerate(cset) if j != idx_in_set && !s.isouter]
                # ck.m_flow
                outerfvs = [get_current_var(namespace, s, fv) for s in cset if s.isouter]
                outersvs = [get_current_var(namespace, s, sv) for s in cset if s.isouter]

                sub[ex] = term(instream_rt, Val(length(innerfvs)), Val(length(outerfvs)), innerfvs..., innersvs..., outerfvs..., outersvs...)
            end
        end
    end

    # additional equations
    additional_eqs = Equation[]
    csets = filter(cset->any(e->e.sys.namespace === namespace, cset.set), csets)
    @show csets
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
            push!(additional_eqs, @show states(cset[1].sys.sys, sv) ~ states(cset[2].sys.sys, sv))
        elseif n_inners == 0 && n_outers == 2
            # we don't expand `instream` in this case.
            v1 = states(cset[1].sys.sys, sv)
            v2 = states(cset[2].sys.sys, sv)
            push!(additional_eqs, v1 ~ instream(v2))
            push!(additional_eqs, v2 ~ instream(v1))
        else
            sq = 0
            s_inners = (s for s in cset if !s.isouter)
            s_outers = (s for s in cset if s.isouter)
            for (q, oscq) in enumerate(s_outers)
                sq += sum(s->max(-states(s, fv), 0), s_inners)
                for (k, s) in enumerate(s_outers); k == q && continue
                    f = states(s.sys.sys, fv)
                    sq += max(f, 0)
                end

                num = 0
                den = 0
                for s in s_inners
                    f = states(s.sys.sys, fv)
                    tmp = positivemax(-f, sq; tol=tol)
                    den += tmp
                    num += tmp * states(s.sys.sys, sv)
                end
                for (k, s) in enumerate(s_outers); k == q && continue
                    f = states(s.sys.sys, fv)
                    tmp = positivemax(f, sq; tol=tol)
                    den += tmp
                    num += tmp * instream(states(s.sys.sys, sv))
                end
                push!(additional_eqs, states(oscq.sys.sys, sv) ~ num / den)
            end
        end
    end
    @show additional_eqs

    display(instream_exprs)
    display(sub)
    @set! sys.systems = []
    @set! sys.eqs = [get_eqs(sys); eqs; substitute(instream_eqs, sub); additional_eqs]
    sys
end

function get_current_var(namespace, cele, sv)
    states(renamespace(unnamespace(namespace, cele.sys.namespace), cele.sys.sys), sv)
end

function get_cset_sv(namespace, ex, csets)
    ns_sv = only(arguments(ex))
    full_name_sv = renamespace(namespace, ns_sv)

    cidx = -1
    idx_in_set = -1
    sv = ns_sv
    for (i, c) in enumerate(csets)
        for (j, v) in enumerate(c.set)
            if isequal(namespaced_var(v), full_name_sv)
                cidx = i
                idx_in_set = j
                sv = v.v
            end
        end
    end
    cidx < 0 && error("$ns_sv is not a variable inside stream connectors")
    cset = csets[cidx].set
    if namespace != first(cset).sys.namespace
        cset = map(c->@set(c.isouter = false), cset)
    end
    cset, idx_in_set, sv
end

@generated function _instream_split(::Val{inner_n}, ::Val{outer_n}, vars::NTuple{N,Any}) where {inner_n, outer_n, N}
    #instream_rt(innerfvs..., innersvs..., outerfvs..., outersvs...)
    ret = Expr(:tuple)
    # mj.c.m_flow
    inner_f = :(Base.@ntuple $inner_n i -> vars[i])
    offset = inner_n
    inner_s = :(Base.@ntuple $inner_n i -> vars[$offset+i])
    offset += inner_n
    # ck.m_flow
    outer_f = :(Base.@ntuple $outer_n i -> vars[$offset+i])
    offset += outer_n
    outer_s = :(Base.@ntuple $outer_n i -> vars[$offset+i])
    Expr(:tuple, inner_f, inner_s, outer_f, outer_s)
end

function instream_rt(ins::Val{inner_n}, outs::Val{outer_n}, vars::Vararg{Any,N}) where {inner_n, outer_n, N}
    @assert N == 2*(inner_n + outer_n)

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

#=
function expand_instream2!(additional_eqs, csets::AbstractVector{<:ConnectionSet}, ins, sys::AbstractSystem, namespace=nothing, prevnamespace=nothing; debug=false, tol=1e-8)
    subsys = get_systems(sys)
    # no connectors if there are no subsystems
    isempty(subsys) && return
    # post order traversal
    for s in subsys
        expand_instream2!(additional_eqs, csets, ins, s, renamespace(namespace, nameof(s)), namespace; debug, tol)
    end

    instream_eqs, instream_exprs = ins[namespace]
    sys = flatten(sys)

    sub = Dict()
    dels = Int[]
    for (k, ex) in enumerate(instream_exprs)
        ex_n = namespace_expr(ex, sys, namespace)
        ns_sv = only(arguments(ex))
        full_name_sv = renamespace(namespace, ns_sv)
        @show full_name_sv

        cidx = -1
        idx_in_set = -1
        sv = ns_sv
        for (i, c) in enumerate(csets)
            for (j, v) in enumerate(c.set)
                if isequal(namespaced_var(v), full_name_sv)
                    cidx = i
                    idx_in_set = j
                    sv = v.v
                end
            end
        end
        #cidx < 0 && error("$ns_sv is not a variable inside stream connectors")
        if cidx > 0
            cset = csets[cidx].set

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
            #@assert all(s->first(cset).sys.namespace == s.sys.namespace, cset)
            if first(cset).sys.namespace != prevnamespace
                cset = map(c->@set(c.isouter = false), cset)
            end
            @show prevnamespace
            @show n_inners, n_outers
            @show nameof(sys) ConnectionSet(cset)
        else
            n_inners = 1
            n_outers = 0
        end
        if n_inners == 1 && n_outers == 0
            @show namespace, ex, sv
            sub[ex_n] = renamespace(renamespace(namespace, nameof(sys)), sv)
        elseif n_inners == 2 && n_outers == 0
            other = idx_in_set == 1 ? 2 : 1
            sub[ex_n] = states(cset[other], sv)
        elseif n_inners == 1 && n_outers == 1
            if !cset[idx_in_set].isouter
                other = idx_in_set == 1 ? 2 : 1
                outerstream = instream(states(cset[other].sys.sys, sv))
                @show outerstream
                ns = nameof(sys)
                ex = namespace_expr(ex, sys, ns)
                # TODO: write a mapping from eqs to exprs
                push!(dels, k)
                eq = substitute(namespace_equation(instream_eqs[k], sys, ns), Dict(ex=>outerstream))
                previnstream_eqs, previnstream_exprs = ins[prevnamespace]
                push!(previnstream_eqs, eq)
                push!(previnstream_exprs, outerstream)
                #outerstream = states(cset[other], sv)
                #substitute()
            end
        else
            if !cset[idx_in_set].isouter
                fv = flowvar(first(connectors))
                # mj.c.m_flow
                innerfvs = [unwrap(states(s, fv)) for (j, s) in enumerate(cset) if j != idx_in_set && !s.isouter]
                innersvs = [unwrap(states(s, sv)) for (j, s) in enumerate(cset) if j != idx_in_set && !s.isouter]
                # ck.m_flow
                outerfvs = [unwrap(states(s, fv)) for s in cset if s.isouter]
                outersvs = [unwrap(states(s, sv)) for s in cset if s.isouter]

                sub[ex_n] = term(instream_rt, Val(length(innerfvs)), Val(length(outerfvs)), innerfvs..., innersvs..., outerfvs..., outersvs...)
            end
        end
    end
    instream_eqs = deleteat!(copy(instream_eqs), dels)
    instream_eqs = map(eq->substitute(namespace_equation(eq, sys, namespace), sub), instream_eqs)

    # additional equations
    csets = filter(cset->any(e->e.sys.namespace === namespace, cset.set), csets)
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
            push!(additional_eqs, states(cset[1], sv) ~ states(cset[2], sv))
        elseif n_inners == 0 && n_outers == 2
            # we don't expand `instream` in this case.
            v1 = states(cset[1], sv)
            v2 = states(cset[2], sv)
            push!(additional_eqs, v1 ~ instream(v2))
            push!(additional_eqs, v2 ~ instream(v1))
        else
            sq = 0
            s_inners = (s for s in cset if !s.isouter)
            s_outers = (s for s in cset if s.isouter)
            for (q, oscq) in enumerate(s_outers)
                sq += sum(s->max(-states(s, fv), 0), s_inners)
                for (k, s) in enumerate(s_outers); k == q && continue
                    f = states(s, fv)
                    sq += max(f, 0)
                end

                num = 0
                den = 0
                for s in s_inners
                    f = states(s, fv)
                    tmp = positivemax(-f, sq; tol=tol)
                    den += tmp
                    num += tmp * states(s, sv)
                end
                for (k, s) in enumerate(s_outers); k == q && continue
                    f = states(s, fv)
                    tmp = positivemax(f, sq; tol=tol)
                    den += tmp
                    num += tmp * instream(states(s, sv))
                end
                push!(additional_eqs, states(oscq, sv) ~ num / den)
            end
        end
    end

    if debug
        println("===========BEGIN=============")
        println("Expanded equations:")
        for eq in instream_eqs
            print_with_indent(4, eq)
        end
        if !isempty(additional_eqs)
            println("Additional equations:")
            for eq in additional_eqs
                print_with_indent(4, eq)
            end
        end
        println("============END==============")
    end
    append!(additional_eqs, instream_eqs)
    return
end
=#

function expand_instream(instream_eqs, instream_exprs, connects; debug=false, tol)
    sub = Dict()
    for ex in instream_exprs
        var = only(arguments(ex))
        connector_name, streamvar_name = split_sys_var(var)

        # find the connect
        cidx = findfirst(c->connector_name in c, connects)
        cidx === nothing && error("$var is not a variable inside stream connectors")
        connect = connects[cidx]
        connectors = Iterators.flatten((connect.inners, connect.outers))
        # stream variable
        sv = getproperty(first(connectors), streamvar_name; namespace=false)
        inner_sc = something(connect.inners, EMPTY_VEC)
        outer_sc = something(connect.outers, EMPTY_VEC)

        n_outers = length(outer_sc)
        n_inners = length(inner_sc)
        outer_names = (nameof(s) for s in outer_sc)
        inner_names = (nameof(s) for s in inner_sc)
        if debug
            println("Expanding: $ex")
            isempty(inner_names) || println("Inner connectors: $(collect(inner_names))")
            isempty(outer_names) || println("Outer connectors: $(collect(outer_names))")
        end

        # expand `instream`s
        # https://specification.modelica.org/v3.4/Ch15.html
        # Based on the above requirements, the following implementation is
        # recommended:
        if n_inners == 1 && n_outers == 0
            connector_name === only(inner_names) || error("$var is not in any stream connector of $(nameof(ogsys))")
            sub[ex_n] = var
        elseif n_inners == 2 && n_outers == 0
            connector_name in inner_names || error("$var is not in any stream connector of $(nameof(ogsys))")
            idx = findfirst(c->nameof(c) === connector_name, inner_sc)
            other = idx == 1 ? 2 : 1
            sub[ex_n] = states(inner_sc[other], sv)
        elseif n_inners == 1 && n_outers == 1
            isinner = connector_name === only(inner_names)
            isouter = connector_name === only(outer_names)
            (isinner || isouter) || error("$var is not in any stream connector of $(nameof(ogsys))")
            if isinner
                outerstream = states(only(outer_sc), sv) # c_1.h_outflow
                sub[ex_n] = outerstream
            end
        else
            fv = flowvar(first(connectors))
            i = findfirst(c->nameof(c) === connector_name, inner_sc)
            if i !== nothing
                # mj.c.m_flow
                innerfvs = [unwrap(states(s, fv)) for (j, s) in enumerate(inner_sc) if j != i]
                innersvs = [unwrap(states(s, sv)) for (j, s) in enumerate(inner_sc) if j != i]
                # ck.m_flow
                outerfvs = [unwrap(states(s, fv)) for s in outer_sc]
                outersvs = [instream(states(s, sv)) for s in outer_sc]

                sub[ex_n] = term(instream_rt, Val(length(innerfvs)), Val(length(outerfvs)), innerfvs..., innersvs..., outerfvs..., outersvs...)
                #=
                si = isempty(outer_sc) ? 0 : sum(s->max(states(s, fv), 0), outer_sc)
                for j in 1:n_inners; j == i && continue
                    f = states(inner_sc[j], fv)
                    si += max(-f, 0)
                end

                num = 0
                den = 0
                for j in 1:n_inners; j == i && continue
                    f = states(inner_sc[j], fv)
                    tmp = positivemax(-f, si; tol=tol)
                    den += tmp
                    num += tmp * states(inner_sc[j], sv)
                end
                for k in 1:n_outers
                    f = states(outer_sc[k], fv)
                    tmp = positivemax(f, si; tol=tol)
                    den += tmp
                    num += tmp * instream(states(outer_sc[k], sv))
                end
                sub[ex_n] = mydiv(num, den)
                =#
            end
        end
    end

    # additional equations
    additional_eqs = Equation[]
    for c in connects
        outer_sc = something(c.outers, EMPTY_VEC)
        isempty(outer_sc) && continue
        inner_sc = c.inners
        n_outers = length(outer_sc)
        n_inners = length(inner_sc)
        connector_representative = first(outer_sc)
        fv = flowvar(connector_representative)
        for sv in get_states(connector_representative)
            vtype = get_connection_type(sv)
            vtype === Stream || continue
            if n_inners == 1 && n_outers == 1
                innerstream = states(only(inner_sc), sv)
                outerstream = states(only(outer_sc), sv)
                push!(additional_eqs, outerstream ~ innerstream)
            elseif n_inners == 0 && n_outers == 2
                # we don't expand `instream` in this case.
                v1 = states(outer_sc[1], sv)
                v2 = states(outer_sc[2], sv)
                push!(additional_eqs, v1 ~ instream(v2))
                push!(additional_eqs, v2 ~ instream(v1))
            else
                sq = 0
                for q in 1:n_outers
                    sq += sum(s->max(-states(s, fv), 0), inner_sc)
                    for k in 1:n_outers; k == q && continue
                        f = states(outer_sc[k], fv)
                        sq += max(f, 0)
                    end

                    num = 0
                    den = 0
                    for j in 1:n_inners
                        f = states(inner_sc[j], fv)
                        tmp = positivemax(-f, sq; tol=tol)
                        den += tmp
                        num += tmp * states(inner_sc[j], sv)
                    end
                    for k in 1:n_outers; k == q && continue
                        f = states(outer_sc[k], fv)
                        tmp = positivemax(f, sq; tol=tol)
                        den += tmp
                        num += tmp * instream(states(outer_sc[k], sv))
                    end
                    push!(additional_eqs, states(outer_sc[q], sv) ~ num / den)
                end
            end
        end
    end

    instream_eqs = map(Base.Fix2(substitute, sub), instream_eqs)
    if debug
        println("===========BEGIN=============")
        println("Expanded equations:")
        for eq in instream_eqs
            print_with_indent(4, eq)
        end
        if !isempty(additional_eqs)
            println("Additional equations:")
            for eq in additional_eqs
                print_with_indent(4, eq)
            end
        end
        println("============END==============")
    end
    return
end
