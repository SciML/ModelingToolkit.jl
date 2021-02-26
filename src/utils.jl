function make_operation(@nospecialize(op), args)
    if op === (*)
        args = filter(!_isone, args)
        if isempty(args)
            return 1
        end
    elseif op === (+)
        args = filter(!_iszero, args)
        if isempty(args)
            return 0
        end
    end
    return op(args...)
end

function detime_dvs(op)
    if !istree(op)
        op
    elseif operation(op) isa Sym
        Sym{Real}(nameof(operation(op)))
    else
        similarterm(op, operation(op),detime_dvs.(arguments(op)))
    end
end

function retime_dvs(op::Sym,dvs,iv)
    Sym{FnType{Tuple{symtype(iv)}, Real}}(nameof(op))(iv)
end

function retime_dvs(op, dvs, iv)
    istree(op) ?
        similarterm(op, operation(op), retime_dvs.(arguments(op),(dvs,),(iv,))) :
        op
end

modified_states!(mstates, e::Equation, statelist=nothing) = get_variables!(mstates, e.lhs, statelist)

# variable substitution
# Piracy but mild
"""
    substitute(expr, s::Pair)
    substitute(expr, s::Dict)
    substitute(expr, s::Vector)

Performs the substitution on `expr` according to rule(s) `s`.

# Examples
```julia
julia> @parameters t
(t,)

julia> @variables x y z(t)
(x, y, z(t))

julia> ex = x + y + sin(z)
(x + y) + sin(z(t))

julia> substitute(ex, Dict([x => z, sin(z) => z^2]))
(z(t) + y) + (z(t) ^ 2)
```
"""
# cannot use Union{Pair, Vector, Dict} -- leads to ambiguity
substitute(expr::Num, s::Pair; kw...) = Num(substituter(s)(value(expr); kw...))
substitute(expr::Num, s::Vector; kw...) = Num(substituter(s)(value(expr); kw...))
substitute(expr::Num, s::Dict; kw...) = Num(substituter(s)(value(expr); kw...))
# TODO: move this to SymbolicUtils
substitute(expr, s::Pair; kw...) = substituter([s[1] => s[2]])(expr; kw...)
substitute(expr, s::Vector; kw...) = substituter(s)(expr; kw...)

substituter(pair::Pair) = substituter((pair,))
function substituter(pairs)
    dict = Dict(value(k) => value(v)  for (k, v) in pairs)
    (expr; kw...) -> SymbolicUtils.substitute(expr, dict; kw...)
end

macro showarr(x)
    n = string(x)
    quote
        y = $(esc(x))
        println($n, " = ", summary(y))
        Base.print_array(stdout, y)
        println()
        y
    end
end

@deprecate substitute_expr!(expr,s) substitute(expr,s)

function states_to_sym(states::Set)
    function _states_to_sym(O)
        if O isa Equation
            Expr(:(=), _states_to_sym(O.lhs), _states_to_sym(O.rhs))
        elseif istree(O)
            op = operation(O)
            args = arguments(O)
            if op isa Sym
                O in states && return tosymbol(O)
                # dependent variables
                return build_expr(:call, Any[nameof(op); _states_to_sym.(args)])
            else
                canonical, O = canonicalexpr(O)
                return canonical ? O : build_expr(:call, Any[op; _states_to_sym.(args)])
            end
        elseif O isa Num
            return _states_to_sym(value(O))
        else
            return toexpr(O)
        end
    end
end
states_to_sym(states) = states_to_sym(Set(states))

"""
    toparam(s::Sym) -> Sym{<:Parameter}

Maps the variable to a paramter.
"""
toparam(s::Sym) = Sym{Parameter{symtype(s)}}(s.name)
toparam(s::Sym{<:Parameter}) = s

"""
    tovar(s::Sym) -> Sym{Real}
    tovar(s::Sym{<:Parameter}) -> Sym{Real}

Maps the variable to a state.
"""
tovar(s::Sym{<:Parameter}) = Sym{symtype(s)}(s.name)
tovar(s::Sym) = s

function lower_mapnames(umap::AbstractArray{T}) where {T<:Pair}
    T[value(k) => value(v) for (k, v) in umap]
end
function lower_mapnames(umap::AbstractArray{T},name) where {T<:Pair}
    T[lower_varname(value(k), name) => value(v) for (k, v) in umap]
end
function lower_mapnames(umap::NTuple{N,T}) where {N,T<:Pair}
    ntuple(i->value(umap[i][1]) => value(umap[i][2]),N)
end
function lower_mapnames(umap::NTuple{N,T},name) where {N,T<:Pair}
    ntuple(i->lower_varname(value(umap[i][1]), name) => value(umap[i][2]),N)
end

lower_mapnames(umap::AbstractArray) = umap # Ambiguity
lower_mapnames(umap::AbstractArray,name) = umap
lower_mapnames(umap::Tuple) = umap
lower_mapnames(umap::Tuple, name) = umap
