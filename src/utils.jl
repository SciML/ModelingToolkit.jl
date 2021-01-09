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

function build_expr(head::Symbol, args)
    ex = Expr(head)
    append!(ex.args, args)
    ex
end

# used in parsing
isblock(x) = length(x) == 1 && x[1] isa Expr && x[1].head == :block
function flatten_expr!(x)
    isb = isblock(x)
    if isb
        x = MacroTools.striplines(x[1])
        filter!(z->z isa Symbol || z.head != :line, x.args)
        x = (x.args...,)
    end
    x
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

is_derivative(O::Term) = isa(operation(O), Differential)
is_derivative(::Any) = false

"""
    get_variables(O) -> Vector{Union{Sym, Term}}

Returns the variables in the expression. Note that the returned variables are
not wrapped in the `Num` type.

# Examples
```julia
julia> @parameters t
(t,)

julia> @variables x y z(t)
(x, y, z(t))

julia> ex = x + y + sin(z)
(x + y) + sin(z(t))

julia> ModelingToolkit.get_variables(ex)
3-element Vector{Any}:
 x
 y
 z(t)
```
"""
get_variables(e::Num, varlist=nothing) = get_variables(value(e), varlist)
get_variables!(vars, e, varlist=nothing) = vars

is_singleton(e::Term) = operation(e) isa Sym
is_singleton(e::Sym) = true
is_singleton(e) = false

get_variables!(vars, e::Number, varlist=nothing) = vars

function get_variables!(vars, e::Symbolic, varlist=nothing)
    if is_singleton(e)
        if isnothing(varlist) || any(isequal(e), varlist)
            push!(vars, e)
        end
    else
        foreach(x -> get_variables!(vars, x, varlist), arguments(e))
    end
    return (vars isa AbstractVector) ? unique!(vars) : vars
end

function get_variables!(vars, e::Equation, varlist=nothing)
  get_variables!(vars, e.rhs, varlist)
end

get_variables(e, varlist=nothing) = get_variables!([], e, varlist)

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
substitute(expr::Num, s::Union{Pair, Vector, Dict}; kw...) = Num(substituter(s)(value(expr); kw...))
# TODO: move this to SymbolicUtils
substitute(expr, s::Pair; kw...) = substituter([s[1] => s[2]])(expr; kw...)
substitute(expr, s::Vector; kw...) = substituter(s)(expr; kw...)

substituter(pair::Pair) = substituter((pair,))
function substituter(pairs)
    dict = Dict(to_symbolic(k) => to_symbolic(v)  for (k, v) in pairs)
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
            if isa(operation(O), Sym)
                O in states && return tosymbol(O)
                # dependent variables
                return build_expr(:call, Any[operation(O).name; _states_to_sym.(arguments(O))])
            else
                return build_expr(:call, Any[operation(O); _states_to_sym.(arguments(O))])
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

Base.Symbol(x::Union{Num,Symbolic}) = tosymbol(x)
tosymbol(x; kwargs...) = x
tosymbol(x::Sym; kwargs...) = nameof(x)
tosymbol(t::Num; kwargs...) = tosymbol(value(t); kwargs...)

"""
    tosymbol(x::Union{Num,Symbolic}; states=nothing, escape=true) -> Symbol

Convert `x` to a symbol. `states` are the states of a system, and `escape`
means if the target has escapes like `val"y⦗t⦘"`. If `escape` then it will only
output `y` instead of `y⦗t⦘`.

# Examples
```julia
julia> @parameters t; @variables z(t)
(z(t),)

julia> ModelingToolkit.tosymbol(z)
Symbol("z⦗t⦘")
```
"""
function tosymbol(t::Term; states=nothing, escape=true)
    if operation(t) isa Sym
        if states !== nothing && !(any(isequal(t), states))
            return nameof(operation(t))
        end
        op = nameof(operation(t))
        args = arguments(t)
    elseif operation(t) isa Differential
        term = diff2term(t)
        op = Symbol(operation(term))
        args = arguments(term)
    else
        @goto err
    end

    return escape ? Symbol(op, "⦗", join(args, ", "), "⦘") : op
    @label err
    error("Cannot convert $t to a symbol")
end

"""
    makesym(x::Union{Num,Symbolic}, kwargs...) -> Sym

`makesym` takes the same arguments as [`tosymbol`](@ref), but it converts a
`Term` in the form of `x(t)` to a `Sym` in the form of `x⦗t⦘`.

# Examples
```julia
julia> @parameters t; @variables x(t)
(x(t),)

julia> ModelingToolkit.makesym(x)
x⦗t⦘
```
"""
makesym(t::Symbolic; kwargs...) = Sym{symtype(t)}(tosymbol(t; kwargs...))
makesym(t::Num; kwargs...) = makesym(value(t); kwargs...)

function lower_varname(var::Symbolic, idv, order)
    order == 0 && return var
    name = string(nameof(operation(var)))
    underscore = 'ˍ'
    idx = findlast(underscore, name)
    append = string(idv)^order
    if idx === nothing
        newname = Symbol(name, underscore, append)
    else
        nidx = nextind(name, idx)
        newname = Symbol(name[1:idx], name[nidx:end], append)
    end
    return Sym{symtype(operation(var))}(newname)(arguments(var)[1])
end

function lower_varname(t::Symbolic, iv)
    var, order = var_from_nested_derivative(t)
    lower_varname(var, iv, order)
end
lower_varname(t::Sym, iv) = t

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

lower_mapnames(umap::AbstractArray{<:Number}) = umap # Ambiguity
lower_mapnames(umap::AbstractArray{<:Number},name) = umap
lower_mapnames(umap::Tuple) = umap
lower_mapnames(umap::Tuple, name) = umap

function flatten_differential(O::Term)
    @assert is_derivative(O) "invalid differential: $O"
    is_derivative(arguments(O)[1]) || return (arguments(O)[1], operation(O).x, 1)
    (x, t, order) = flatten_differential(arguments(O)[1])
    isequal(t, operation(O).x) || throw(ArgumentError("non-matching differentials on lhs: $t, $(operation(O).x)"))
    return (x, t, order + 1)
end

"""
    diff2term(x::Term) -> Term
    diff2term(x) -> x

Convert a differential variable to a `Term`. Note that it only takes a `Term`
not a `Num`.
```julia
julia> ModelingToolkit.diff2term(ModelingToolkit.value(D(D(x))))
xˍtt(t)
```
"""
function diff2term(O)
    istree(O) || return O
    if is_derivative(O)
        (x, t, order) = flatten_differential(O)
        return lower_varname(x, t, order)
    end
    return Term{Real}(operation(O), diff2term.(arguments(O)))
end
