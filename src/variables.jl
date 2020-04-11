const IndexMap = Dict{Char,Char}(
            '0' => '₀',
            '1' => '₁',
            '2' => '₂',
            '3' => '₃',
            '4' => '₄',
            '5' => '₅',
            '6' => '₆',
            '7' => '₇',
            '8' => '₈',
            '9' => '₉')
function map_subscripts(indices)
    str = string(indices)
    join(IndexMap[c] for c in str)
end

"""
$(TYPEDEF)

A named variable which represents a numerical value.

# Fields
$(FIELDS)
"""
struct Variable{T} <: Function
    """The variable's unique name."""
    name::Symbol
    Variable(name) = new{Number}(name)
    Variable{T}(name) where T = new{T}(name)
end

getindex(x::Variable{<:Number},idxs...) = error("Variable $(x.name) is a scalar quantity. No indexing allowed")
getindex(x::Variable{T},idxs...) where T = Variable{eltype(T)}(varSymbol("$(x.name)$(join(map_subscripts.(idxs), "ˏ"))"))

vartype(::Variable{T}) where T = T
(x::Variable)(args...) = Operation(x, collect(Expression, args))

Base.isequal(x::Variable, y::Variable) = x.name == y.name
Base.print(io::IO, x::Variable) = show(io, x)
Base.show(io::IO, x::Variable) = print(io, x.name)
function Base.show(io::IO, ::MIME"text/plain", x::Variable)
    print(io, x.name)
end


"""
$(TYPEDEF)

An expression which wraps a constant numerical value.

The value of the constant can be extracted with [`Base.get`](@ref).
"""
struct Constant <: Expression
    """The constant's numerical value"""
    value::Number
end

"""
$(TYPEDSIGNATURES)

Get the value of a [`ModelingToolkit.Constant`](@ref).
"""
Base.get(c::Constant) = c.value

Base.iszero(ex::Expression) = isa(ex, Constant) && iszero(ex.value)
Base.isone(ex::Expression)  = isa(ex, Constant) && isone(ex.value)

# Variables use isequal for equality since == is an Operation
Base.isequal(c::Constant, n::Number) = c.value == n
Base.isequal(n::Number, c::Constant) = c.value == n
Base.isequal(a::Constant, b::Constant) = a.value == b.value

Base.convert(::Type{Expr}, c::Constant) = c.value


# Build variables more easily
function _parse_vars(macroname, type, x)
    ex = Expr(:block)
    var_names = Symbol[]
    # if parsing things in the form of
    # begin
    #     x
    #     y
    #     z
    # end
    x = x isa Tuple && first(x) isa Expr && first(x).head == :tuple ? first(x).args : x # tuple handling
    x = flatten_expr!(x)
    for _var in x
        iscall = isa(_var, Expr) && _var.head == :call
        isarray = isa(_var, Expr) && _var.head == :ref
        issym  = _var isa Symbol
        @assert iscall || isarray || issym "@$macroname expects a tuple of expressions or an expression of a tuple (`@$macroname x y z(t) v[1:3] w[1:2,1:4]` or `@$macroname x, y, z(t) v[1:3] w[1:2,1:4]`)"

        if iscall
            var_name, expr = _construct_vars(_var.args[1], type, _var.args[2:end])
        else
            var_name, expr = _construct_vars(_var, type, nothing)
        end
        push!(var_names, var_name)
        push!(ex.args, expr)
    end
    push!(ex.args, build_expr(:tuple, var_names))
    return ex
end

function _construct_vars(_var, type, call_args)
    issym  = _var isa Symbol
    isarray = isa(_var, Expr) && _var.head == :ref
    if isarray
        var_name = _var.args[1]
        indices = _var.args[2:end]
        expr = _construct_array_vars(var_name, type, call_args, indices...)
    else
        # Implicit 0-args call
        var_name = _var
        expr = _construct_var(var_name, type, call_args)
    end
    var_name, :($var_name = $expr)
end

function _construct_var(var_name, type, call_args)
    if call_args === nothing
        :(Variable{$type}($(Meta.quot(var_name)))())
    elseif !isempty(call_args) && call_args[end] == :..
        :(Variable{$type}($(Meta.quot(var_name))))
    else
        :(Variable{$type}($(Meta.quot(var_name)))($(call_args...)))
    end
end

function _construct_var(var_name, type, call_args, ind)
    if call_args === nothing
        :(Variable{$type}($(Meta.quot(var_name)), $ind...)())
    elseif !isempty(call_args) && call_args[end] == :..
        :(Variable{$type}($(Meta.quot(var_name)), $ind...))
    else
        :(Variable{$type}($(Meta.quot(var_name)), $ind...)($(call_args...)))
    end
end


function _construct_array_vars(var_name, type, call_args, indices...)
    T = :(SizedArray{Tuple{length.($(indices))...},$type})
    @show T
    _construct_var(var_name, T, call_args)
end


"""
$(SIGNATURES)

Define one or more unknown variables.
"""
macro variables(xs...)
    esc(_parse_vars(:variables, Number, xs))
end

function rename(x::Variable,name::Symbol)
    Variable{vartype(x)}(name)
end

TreeViews.hastreeview(x::Variable) = true
function TreeViews.treelabel(io::IO,x::Variable,
                             mime::MIME"text/plain" = MIME"text/plain"())
  show(io,mime,Text(x.name))
end
