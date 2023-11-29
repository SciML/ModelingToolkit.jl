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
    if istree(s) && operation(s) === getindex
        s = arguments(s)[1]
    end
    getmetadata(s, VariableConnectType, Equality)
end

macro connector(expr)
    esc(component_post_processing(expr, true))
end

abstract type AbstractConnectorType end
struct StreamConnector <: AbstractConnectorType end
struct RegularConnector <: AbstractConnectorType end
struct DomainConnector <: AbstractConnectorType end

function connector_type(sys::AbstractSystem)
    sts = get_states(sys)
    n_stream = 0
    n_flow = 0
    n_regular = 0 # state that is not input, output, stream, or flow.
    for s in sts
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
    if n_flow == 1 && length(sts) == 1
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

instream(a) = term(instream, unwrap(a), type = symtype(a))
SymbolicUtils.promote_symtype(::typeof(instream), _) = Real

isconnector(s::AbstractSystem) = has_connector_type(s) && get_connector_type(s) !== nothing

function flowvar(sys::AbstractSystem)
    sts = get_states(sys)
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
            idx = findfirst(isequal('₊'), s)
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
end
Base.nameof(l::ConnectionElement) = renamespace(nameof(l.sys), getname(l.v))
function Base.hash(l::ConnectionElement, salt::UInt)
    hash(nameof(l.sys)) ⊻ hash(l.v) ⊻ hash(l.isouter) ⊻ salt
end
Base.isequal(l1::ConnectionElement, l2::ConnectionElement) = l1 == l2
function Base.:(==)(l1::ConnectionElement, l2::ConnectionElement)
    nameof(l1.sys) == nameof(l2.sys) && isequal(l1.v, l2.v) && l1.isouter == l2.isouter
end
namespaced_var(l::ConnectionElement) = states(l, l.v)
states(l::ConnectionElement, v) = states(copy(l.sys), v)

struct ConnectionSet
    set::Vector{ConnectionElement} # namespace.sys, var, isouter
end
Base.copy(c::ConnectionSet) = ConnectionSet(copy(c.set))

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

"Return true if the system is a 3D multibody frame, otherwise return false."
function isframe(sys)
    (has_metadata(sys) && (md = get_metadata(sys)) isa Dict) || return false
    get(md, :frame, false)
end

"Return orientation object of a multibody frame."
function ori(sys)
    @assert has_metadata(sys)
    md = get_metadata(sys)
    if md isa Dict && (O = get(md, :orientation, nothing)) !== nothing
        return O
    else
        error("System $(sys.name) does not have an orientation object.")
    end
end

function connection2set!(connectionsets, namespace, ss, isouter)
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
        dv = only(states(domain_ss))
        for (i, s) in enumerate(ss)
            sts = states(s)
            io = isouter(s)
            for (j, v) in enumerate(sts)
                vtype = get_connection_type(v)
                (vtype === Flow && isequal(v, dv)) || continue
                push!(cset, T(LazyNamespace(namespace, domain_ss), dv, false))
                push!(cset, T(LazyNamespace(namespace, s), v, io))
            end
        end
        @assert length(cset) > 0
        push!(connectionsets, ConnectionSet(cset))
        return connectionsets
    end
    s1 = first(ss)
    sts1v = states(s1)
    if isframe(s1) # Multibody
        O = ori(s1)
        orientation_vars = Symbolics.unwrap.(collect(vec(O.R)))
        sts1v = [sts1v; orientation_vars]
    end
    sts1 = Set(sts1v)
    num_statevars = length(sts1)
    csets = [T[] for _ in 1:num_statevars] # Add 9 orientation variables if connection is between multibody frames
    for (i, s) in enumerate(ss)
        sts = states(s)
        if isframe(s) # Multibody
            O = ori(s)
            orientation_vars = Symbolics.unwrap.(vec(O.R))
            sts = [sts; orientation_vars]
        end
        i != 1 && ((num_statevars == length(sts) && all(Base.Fix2(in, sts1), sts)) ||
         connection_error(ss))
        io = isouter(s)
        for (j, v) in enumerate(sts)
            push!(csets[j], T(LazyNamespace(namespace, s), v, io))
        end
    end
    for cset in csets
        v = first(cset).v
        vtype = get_connection_type(v)
        if domain_ss !== nothing && vtype === Flow &&
           (dv = only(states(domain_ss)); isequal(v, dv))
            push!(cset, T(LazyNamespace(namespace, domain_ss), dv, false))
        end
        for k in 2:length(cset)
            vtype === get_connection_type(cset[k].v) || connection_error(ss)
        end
        push!(connectionsets, ConnectionSet(cset))
    end
end

function generate_connection_set(sys::AbstractSystem, find = nothing, replace = nothing)
    connectionsets = ConnectionSet[]
    domain_csets = ConnectionSet[]
    sys = generate_connection_set!(connectionsets, domain_csets, sys, find, replace)
    csets = merge(connectionsets)
    domain_csets = merge([csets; domain_csets], true)

    sys, (csets, domain_csets)
end

function generate_connection_set!(connectionsets, domain_csets,
        sys::AbstractSystem, find, replace, namespace = nothing)
    subsys = get_systems(sys)

    isouter = generate_isouter(sys)
    eqs′ = get_eqs(sys)
    eqs = Equation[]

    cts = [] # connections
    domain_cts = [] # connections
    extra_states = []
    for eq in eqs′
        lhs = eq.lhs
        rhs = eq.rhs
        if find !== nothing && find(rhs, _getname(namespace))
            neweq, extra_state = replace(rhs, _getname(namespace))
            if extra_state isa AbstractArray
                append!(extra_states, unwrap.(extra_state))
            elseif extra_state !== nothing
                push!(extra_states, extra_state)
            end
            neweq isa AbstractArray ? append!(eqs, neweq) : push!(eqs, neweq)
        else
            if lhs isa Number || lhs isa Symbolic
                push!(eqs, eq) # split connections and equations
            elseif lhs isa Connection && get_systems(lhs) === :domain
                connection2set!(domain_csets, namespace, get_systems(rhs), isouter)
            else
                push!(cts, get_systems(rhs))
            end
        end
    end

    # all connectors are eventually inside connectors.
    T = ConnectionElement
    for s in subsys
        isconnector(s) || continue
        is_domain_connector(s) && continue
        for v in states(s)
            Flow === get_connection_type(v) || continue
            push!(connectionsets, ConnectionSet([T(LazyNamespace(namespace, s), v, false)]))
        end
    end

    for ct in cts
        connection2set!(connectionsets, namespace, ct, isouter)
    end

    # pre order traversal
    if !isempty(extra_states)
        @set! sys.states = [get_states(sys); extra_states]
    end
    @set! sys.systems = map(s -> generate_connection_set!(connectionsets, domain_csets, s,
            find, replace,
            renamespace(namespace, s)),
        subsys)
    @set! sys.eqs = eqs
end

function Base.merge(csets::AbstractVector{<:ConnectionSet}, allouter = false)
    csets, merged = partial_merge(csets, allouter)
    while merged
        csets, merged = partial_merge(csets)
    end
    csets
end

function partial_merge(csets::AbstractVector{<:ConnectionSet}, allouter = false)
    mcsets = ConnectionSet[]
    ele2idx = Dict{ConnectionElement, Int}()
    cacheset = Set{ConnectionElement}()
    merged = false
    for (j, cset) in enumerate(csets)
        if allouter
            cset = ConnectionSet(map(cset.set) do e
                @set! e.isouter = true
            end)
        end
        idx = nothing
        for e in cset.set
            idx = get(ele2idx, e, nothing)
            if idx !== nothing
                merged = true
                break
            end
        end
        if idx === nothing
            push!(mcsets, copy(cset))
            for e in cset.set
                ele2idx[e] = length(mcsets)
            end
        else
            for e in mcsets[idx].set
                push!(cacheset, e)
            end
            for e in cset.set
                push!(cacheset, e)
            end
            empty!(mcsets[idx].set)
            for e in cacheset
                ele2idx[e] = idx
                push!(mcsets[idx].set, e)
            end
            empty!(cacheset)
        end
    end
    mcsets, merged
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
                ns_s_def = Dict(states(m.sys.sys, n) => n for (n, v) in s_def)
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

function expand_connections(sys::AbstractSystem, find = nothing, replace = nothing;
        debug = false, tol = 1e-10)
    sys, (csets, domain_csets) = generate_connection_set(sys, find, replace)
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
    @set! sys.systems = map(s -> expand_instream(csets, s,
            renamespace(namespace, nameof(s)),
            namespace; debug, tol), subsys)
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
            push!(additional_eqs, states(cset[1].sys.sys, sv) ~ states(cset[2].sys.sys, sv))
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
                sq += sum(s -> max(-states(s, fv), 0), s_inners, init = 0)
                for (k, s) in enumerate(s_outers)
                    k == q && continue
                    f = states(s.sys.sys, fv)
                    sq += max(f, 0)
                end

                num = 0
                den = 0
                for s in s_inners
                    f = states(s.sys.sys, fv)
                    tmp = positivemax(-f, sq; tol = tol)
                    den += tmp
                    num += tmp * states(s.sys.sys, sv)
                end
                for (k, s) in enumerate(s_outers)
                    k == q && continue
                    f = states(s.sys.sys, fv)
                    tmp = positivemax(f, sq; tol = tol)
                    den += tmp
                    num += tmp * instream(states(s.sys.sys, sv))
                end
                push!(additional_eqs, states(oscq.sys.sys, sv) ~ num / den)
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
    states(renamespace(unnamespace(namespace, _getname(cele.sys.namespace)), cele.sys.sys),
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
