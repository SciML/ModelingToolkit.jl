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
    sts = states(sys)
    #TODO: check the criteria for stream connectors
    any(s->getmetadata(s, ModelingToolkit.VariableConnectType, nothing) === Stream, sts) ? StreamConnector() : RegularConnector()
end

Base.@kwdef struct Connection
    inners = nothing
    outers = nothing
end

# everything is inner by default until we expand the connections
Connection(syss) = Connection(inners=syss)
get_systems(c::Connection) = c.inners

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

function connect(syss::AbstractSystem...)
    length(syss) >= 2 || error("connect takes at least two systems!")
    length(unique(nameof, syss)) == length(syss) || error("connect takes distinct systems!")
    Equation(Connection(), Connection(syss)) # the RHS are connected systems
end

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

isconnector(s::AbstractSystem) = has_connector_type(s) && get_connector_type(s) !== nothing

function isouterconnector(sys::AbstractSystem; check=true)
    subsys = get_systems(sys)
    outer_connectors = [nameof(s) for s in subsys if isconnector(s)]
    # Note that subconnectors in outer connectors are still outer connectors.
    # Ref: https://specification.modelica.org/v3.4/Ch9.html see 9.1.2
    let outer_connectors=outer_connectors, check=check
        function isouter(sys)::Bool
            s = string(nameof(sys))
            check && (isconnector(sys) || error("$s is not a connector!"))
            idx = findfirst(isequal('₊'), s)
            parent_name = Symbol(idx === nothing ? s : s[1:idx])
            parent_name in outer_connectors
        end
    end
end

print_with_indent(n, x) = println(" " ^ n, x)

function expand_connections(sys::AbstractSystem; debug=false)
    subsys = get_systems(sys)
    isempty(subsys) && return sys

    # post order traversal
    @set! sys.systems = map(s->expand_connections(s, debug=debug), subsys)
    isouter = isouterconnector(sys)

    eqs′ = get_eqs(sys)
    eqs = Equation[]
    cts = [] # connections
    for eq in eqs′
        eq.lhs isa Connection ? push!(cts, get_systems(eq.rhs)) : push!(eqs, eq) # split connections and equations
    end

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

    # Bad things happen when there are more than one intersections
    for c in narg_connects
        @unpack outers, inners = c
        len = length(outers) + length(inners)
        allconnectors = Iterators.flatten((outers, inners))
        dups = find_duplicates(nameof(c) for c in allconnectors)
        length(dups) == 0 || error("$(Connection(syss)) has duplicated connections: $(dups).")
    end

    if debug && !isempty(narg_connects)
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

    if debug && !isempty(narg_connects)
        println("Connection equations:")
        foreach(Base.Fix1(print_with_indent, 4), connection_eqs)
        println("=============END=================")
    end

    @set! sys.eqs = eqs
    return sys
end
