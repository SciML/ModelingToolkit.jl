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
            $isdefined(res, :connector_type) ? $Setfield.@set!(res.connector_type = $connector_type(res)) : res
        end
    end
end

macro connector(expr)
    esc(with_connector_type(expr))
end

function connector_type(sys::AbstractSystem)
    states(sys)
end

promote_connect_rule(::Type{T}, ::Type{S}) where {T, S} = Union{}
promote_connect_rule(::Type{T}, ::Type{T}) where {T} = T
promote_connect_type(t1::Type, t2::Type, ts::Type...) = promote_connect_type(promote_connect_rule(t1, t2), ts...)
@inline function promote_connect_type(::Type{T}, ::Type{S}) where {T,S}
    promote_connect_result(
        T,
        S,
        promote_connect_rule(T,S),
        promote_connect_rule(S,T)
    )
end

promote_connect_result(::Type, ::Type, ::Type{T}, ::Type{Union{}}) where {T} = T
promote_connect_result(::Type, ::Type, ::Type{Union{}}, ::Type{S}) where {S} = S
promote_connect_result(::Type, ::Type, ::Type{T}, ::Type{T}) where {T} = T
function promote_connect_result(::Type{T}, ::Type{S}, ::Type{P1}, ::Type{P2}) where {T,S,P1,P2}
    throw(ArgumentError("connection promotion for $T and $S resulted in $P1 and $P2. " *
                        "Define promotion only in one direction."))
end

throw_connector_promotion(T, S) = throw(ArgumentError("Don't know how to connect systems of type $S and $T"))
promote_connect_result(::Type{T},::Type{S},::Type{Union{}},::Type{Union{}}) where {T,S} = throw_connector_promotion(T,S)

promote_connect_type(::Type{T}, ::Type{T}) where {T} = T
function promote_connect_type(T, S)
    error("Don't know how to connect systems of type $S and $T")
end

Base.@kwdef struct Connection
    inners = nothing
    outers = nothing
end

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
        sys_str = join((string(nameof(s)) * (i <= splitting_idx ? ("::inner") : ("::outers")) for (i, s) in enumerate(syss)), ", ")
        print(io, "<", sys_str, ">")
    end
end

function connect(syss...)
    length(syss) >= 2 || error("connect takes at least two systems!")
    length(unique(nameof, syss)) == length(syss) || error("connect takes distinct systems!")
    Equation(Connection(), Connection(syss)) # the RHS are connected systems
end

isconnector(s::AbstractSystem) = has_connector_type(s) && get_connector_type(s) !== nothing

function isouterconnector(sys::AbstractSystem; check=true)
    subsys = get_systems(sys)
    outer_connectors = [nameof(s) for s in subsys if isconnector(sys)]
    # Note that subconnectors in outer connectors are still outer connectors.
    # Ref: https://specification.modelica.org/v3.4/Ch9.html see 9.1.2
    let outer_connectors=outer_connectors, check=check
        function isouter(sys)::Bool
            s = string(nameof(sys))
            check && (isconnector(sys) || error("$s is not a connector!"))
            idx = findfirst(isequal('₊'), s)
            parent_name = idx === nothing ? s : s[1:idx]
            parent_name in outer_connectors
        end
    end
end

function expand_connections(sys::AbstractSystem; debug=false)
    subsys = get_systems(sys)
    isempty(subsys) && return sys

    # post order traversal
    @set sys.systems = map(s->expand_connections(s, debug=debug), subsys)
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

    if debug
        println("Connections:")
        print_with_indent(x) = println(" " ^ 4, x)
        foreach(print_with_indent, narg_connects)
    end

    for c in narg_connects
        T = promote_connect_type(map(get_connector_type, c.outers)..., map(get_connector_type, c.inners)...)
        ceqs = connect(T, c)
        ceqs isa Equation ? push!(eqs, ceqs) : append!(eqs, ceqs)
    end

    @set! sys.eqs = eqs
    return sys
end
