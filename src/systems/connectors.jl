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
        vtype = getmetadata(s, ModelingToolkit.VariableConnectType, nothing)
        vtype === Stream && (n_stream += 1)
        vtype === Flow   && (n_flow += 1)
    end
    (n_stream > 1 && n_flow > 1) && error("There are multiple flow variables in $(nameof(sys))!")
    n_stream > 1 ? StreamConnector() : RegularConnector()
end

Base.@kwdef struct Connection
    inners = nothing
    outers = nothing
end

# everything is inner by default until we expand the connections
Connection(syss) = Connection(inners=syss)
get_systems(c::Connection) = c.inners
function Base.in(e::Symbol, c::Connection)
    any(k->nameof(k) === e, c.inners) || any(k->nameof(k) === e, c.outers)
end

const EMPTY_VEC = []

function Base.show(io::IO, c::Connection)
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
function connect(syss::AbstractSystem...)
    length(syss) >= 2 || error("connect takes at least two systems!")
    length(unique(nameof, syss)) == length(syss) || error("connect takes distinct systems!")
    Equation(Connection(), Connection(syss)) # the RHS are connected systems
end

# the actual `connect`.
function connect(c::Connection; check=true)
    @unpack inners, outers = c

    flow_eqs = Equation[]
    other_eqs = Equation[]

    cnts = Iterators.flatten((inners, outers))
    fs, ss = Iterators.peel(cnts)
    splitting_idx = length(inners) # anything after the splitting_idx is outer.
    first_sts = get_states(fs)
    first_sts_set = Set(getname.(first_sts))
    for sys in ss
        current_sts = getname.(get_states(sys))
        Set(current_sts) == first_sts_set || error("$(nameof(sys)) ($current_sts) doesn't match the connection type of $(nameof(fs)) ($first_sts).")
    end

    ceqs = Equation[]
    for s in first_sts
        name = getname(s)
        isflow = getmetadata(s, VariableConnectType, Equality) === Flow
        rhs = 0 # only used for flow variables
        fix_val = getproperty(fs, name) # used for equality connections
        for (i, c) in enumerate(cnts)
            isinner = i <= splitting_idx
            # https://specification.modelica.org/v3.4/Ch15.html
            var = getproperty(c, name)
            if isflow
                rhs += isinner ? var : -var
            else
                i == 1 && continue # skip the first iteration
                push!(ceqs, fix_val ~ getproperty(c, name))
            end
        end
        isflow && push!(ceqs, 0 ~ rhs)
    end

    return ceqs
end

instream(a) = term(instream, unwrap(a), type=symtype(a))

isconnector(s::AbstractSystem) = has_connector_type(s) && get_connector_type(s) !== nothing
isstreamconnector(s::AbstractSystem) = isconnector(s) && get_connector_type(s) === Stream

print_with_indent(n, x) = println(" " ^ n, x)

function get_stream_connectors!(sc, sys::AbstractSystem)
    subsys = get_systems(sys)
    isempty(subsys) && return nothing
    for s in subsys; isstreamconnector(s) || continue
        push!(sc, renamespace(sys, s))
    end
    for s in subsys
        get_stream_connectors!(sc, renamespace(sys, s))
    end
    nothing
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

function split_var(var)
    name = string(nameof(var))
    map(Symbol, split(name, '₊'))
end

# inclusive means the first level of `var` is `sys`
function get_sys_var(sys::AbstractSystem, var; inclusive=true)
    lvs = split_var(var)
    if inclusive
        sysn, lvs = Iterator.peel(lvs)
        sysn === nameof(sys) || error("$(nameof(sys)) doesn't have $var!")
    end
    newsys = getproperty(sys, first(lvs))
    for i in 2:length(lvs)-1
        newsys = getproperty(newsys, lvs[i])
    end
    newsys, lvs[end]
end

function split_sys_var(var)
    var_name = string(getname(var))
    sidx = findlast(isequal('₊'), var_name)
    sidx === nothing && error("$var is not a namespaced variable")
    connector_name = Symbol(var_name[1:prevind(var_name, sidx)])
    streamvar_name = Symbol(var_name[nextind(var_name, sidx):end])
    connector_name, streamvar_name
end

function find_connection(connector_name, ogsys, names)
    cs = get_connections(ogsys)
    cs === nothing || for c in cs
        renamespace(names, connector_name) in c && return ogsys, c
    end
    innersys = ogsys
    for (i, n) in enumerate(names)
        innersys = getproperty(innersys, n)
        cs = get_connections(innersys)
        cs === nothing || for c in cs
            nn = @view names[i+1:end]
            renamespace(nn, connector_name) in c && return innersys, c
        end
    end
    error("$connector_name cannot be found in $(nameof(ogsys)) with levels $(names)")
end

function flowvar(sys::AbstractSystem)
    sts = get_states(sys)
    for s in sts
        vtype = getmetadata(s, ModelingToolkit.VariableConnectType, nothing)
        vtype === Flow && return s
    end
    error("There in no flow variable in $(nameof(sys))")
end

function expand_instream(ogsys, sys::AbstractSystem=ogsys, names=[]; debug=false)
    subsys = get_systems(sys)
    isempty(subsys) && return sys

    # post order traversal
    @set! sys.systems = map(subsys) do s
        n = nameof(s)
        expand_instream(ogsys, s, [names; n], debug=debug)
    end
    eqs′ = get_eqs(sys)
    eqs = Equation[]
    instream_eqs = Equation[]
    instream_exprs = Set()
    for eq in eqs′
        if collect_instream!(instream_exprs, eq)
            push!(instream_eqs, eq)
        else
            push!(eqs, eq) # split instreams and equations
        end
    end
    #@show nameof(sys), names, instream_exprs
    isempty(instream_eqs) && return sys

    sub = Dict()
    seen = Set()
    additional_eqs = Equation[]
    for ex in instream_exprs
        var = only(arguments(ex))
        connector_name, streamvar_name = split_sys_var(var)

        outer_sc = []
        inner_sc = []
        # find the connect
        parentsys, connect = find_connection(connector_name, ogsys, names)
        connectors = Iterators.flatten((connect.inners, connect.outers))
        # stream variable
        sv = getproperty(first(connectors), streamvar_name; namespace=false)
        if nameof(parentsys) != nameof(sys)
            # everything is a inner connector w.r.t. `sys`
            for s in connectors
                push!(inner_sc, s)
            end
        else
            for s in connectors
                if connector_name == split_var(nameof(s))[1]
                    push!(inner_sc, s)
                else
                    push!(outer_sc, s)
                end
            end
        end

        n_inners = length(outer_sc)
        n_outers = length(inner_sc)
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
            connector_name === only(inner_names) || error("$stream_name is not in any stream connector of $(nameof(ogsys))")
            sub[ex] = var
        elseif n_inners == 2 && n_outers == 0
            connector_name in inner_names || error("$stream_name is not in any stream connector of $(nameof(ogsys))")
            idx = findfirst(c->nameof(c) === connector_name, inner_sc)
            other = idx == 1 ? 2 : 1
            sub[ex] = states(inner_sc[other], sv)
        elseif n_inners == 1 && n_outers == 1
            isinner = connector_name === only(inner_names)
            isouter = connector_name === only(outer_names)
            (isinner || isouter) || error("$stream_name is not in any stream connector of $(nameof(ogsys))")
            if isinner
                outerinstream = states(only(outer_sc), sv) # c_1.h_outflow
                sub[ex] = outerinstream
            end
            if var in seen
                push!(additional_eqs, outerinstream ~ var)
                push!(seen, var)
            end
        elseif n_inners == 0 && n_outers == 2
            # we don't expand `instream` in this case.
            if var in seen
                v1 = states(outer_sc[1], sv)
                v2 = states(outer_sc[2], sv)
                push!(additional_eqs, v1 ~ instream(v2))
                push!(additional_eqs, v2 ~ instream(v1))
                push!(seen, var)
            end
        else
            fv = flowvar(first(connectors))
            idx = findfirst(c->nameof(c) === connector_name, inner_sc)
            if idx !== nothing
                si = sum(s->max(states(s, fv), 0), outer_sc)
                for j in 1:n_inners; j == i && continue
                    f = states(inner_sc[j], fv)
                    si += max(-f, 0)
                end

                num = 0
                den = 0
                for j in 1:n_inners; j == i && continue
                    f = states(inner_sc[j], fv)
                    tmp = positivemax(-f, si)
                    den += tmp
                    num += tmp * states(inner_sc[j], sv)
                end
                for k in 1:n_outers
                    f = states(outer_sc[k], fv)
                    tmp = positivemax(f, si)
                    den += tmp
                    num += tmp * instream(states(outer_sc[k], sv))
                end
                sub[ex] = num / den
            end

            if var in seen
                for q in 1:n_outers
                    sq += sum(s->max(-states(s, fv), 0), inner_sc)
                    for k in 1:n_outers; k == q && continue
                        f = states(outer_sc[j], fv)
                        si += max(f, 0)
                    end

                    num = 0
                    den = 0
                    for j in 1:n_inners
                        f = states(inner_sc[j], fv)
                        tmp = positivemax(-f, sq)
                        den += tmp
                        num += tmp * states(inner_sc[j], sv)
                    end
                    for k in 1:n_outers; k == q && continue
                        f = states(outer_sc[k], fv)
                        tmp = positivemax(f, sq)
                        den += tmp
                        num += tmp * instream(states(outer_sc[k], sv))
                    end
                    push!(additional_eqs, states(outer_sc[q], sv) ~ num / den)
                end
                push!(seen, var)
            end
        end
    end
    instream_eqs = map(Base.Fix2(substitute, sub), instream_eqs)
    if debug
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
    end
    @set! sys.eqs = [eqs; instream_eqs; additional_eqs]
    return flatten(sys)
end

function expand_connections(sys::AbstractSystem; debug=false)
    sys = collect_connections(sys; debug=debug)
    sys = expand_instream(sys; debug=debug)
    return sys
end

function collect_connections(sys::AbstractSystem; debug=false)
    subsys = get_systems(sys)
    isempty(subsys) && return sys

    # post order traversal
    @set! sys.systems = map(s->collect_connections(s, debug=debug), subsys)

    outer_connectors = Symbol[]
    for s in subsys
        n = nameof(s)
        isconnector(s) && push!(outer_connectors, n)
    end
    isouter = let outer_connectors=outer_connectors
        function isouter(sys)::Bool
            s = string(nameof(sys))
            isconnector(sys) || error("$s is not a connector!")
            idx = findfirst(isequal('₊'), s)
            parent_name = Symbol(idx === nothing ? s : s[1:prevind(s, idx)])
            parent_name in outer_connectors
        end
    end

    eqs′ = get_eqs(sys)
    eqs = Equation[]
    cts = [] # connections
    for eq in eqs′
        eq.lhs isa Connection ? push!(cts, get_systems(eq.rhs)) : push!(eqs, eq) # split connections and equations
    end

    # if there are no connections, we are done
    isempty(cts) && return sys

    sys2idx = Dict{Symbol,Int}() # system (name) to n-th connect statement
    narg_connects = Connection[]
    for (i, syss) in enumerate(cts)
        # find intersecting connections
        exclude = 0 # exclude the intersecting system
        idx = 0     # idx of narg_connects
        for (j, s) in enumerate(syss)
            idx′ = get(sys2idx, nameof(s), nothing)
            idx′ === nothing && continue
            idx = idx′
            exclude = j
        end
        if exclude == 0
            outers = []
            inners = []
            for s in syss
                isouter(s) ? push!(outers, s) : push!(inners, s)
            end
            push!(narg_connects, Connection(outers=outers, inners=inners))
            for s in syss
                sys2idx[nameof(s)] = length(narg_connects)
            end
        else
            # fuse intersecting connections
            for (j, s) in enumerate(syss); j == exclude && continue
                sys2idx[nameof(s)] = idx
                c = narg_connects[idx]
                isouter(s) ? push!(c.outers, s) : push!(c.inners, s)
            end
        end
    end

    isempty(narg_connects) && error("Unreachable reached. Please file an issue.")

    @set! sys.connections = narg_connects

    # Bad things happen when there are more than one intersections
    for c in narg_connects
        @unpack outers, inners = c
        len = length(outers) + length(inners)
        allconnectors = Iterators.flatten((outers, inners))
        dups = find_duplicates(nameof(c) for c in allconnectors)
        length(dups) == 0 || error("$(Connection(syss)) has duplicated connections: $(dups).")
    end

    if debug
        println("============BEGIN================")
        println("Connections for [$(nameof(sys))]:")
        foreach(Base.Fix1(print_with_indent, 4), narg_connects)
    end

    connection_eqs = Equation[]
    for c in narg_connects
        ceqs = connect(c)
        debug && append!(connection_eqs, ceqs)
        append!(eqs, ceqs)
    end

    if debug
        println("Connection equations:")
        foreach(Base.Fix1(print_with_indent, 4), connection_eqs)
        println("=============END=================")
    end

    @set! sys.eqs = eqs
    return sys
end
