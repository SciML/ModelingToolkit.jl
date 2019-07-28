export Variable, @variables, @parameters

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

A named variable which represents a numerical value. The variable's value may
be known (parameters, independent variables) or unknown (dependent variables).

# Fields
$(FIELDS)
"""
struct Variable <: Function
    """The variable's unique name."""
    name::Symbol
    """
    Whether the variable's value is known.
    """
    known::Bool
    Variable(name; known = false) = new(name, known)
end
function Variable(name, indices...; known = false)
    var_name = Symbol("$(name)$(join(map_subscripts.(indices), "̒"))")
    Variable(var_name; known=known)
end

(x::Variable)(args...) = Operation(x, collect(Expression, args))

Base.isequal(x::Variable, y::Variable) = (x.name, x.known) == (y.name, y.known)
Base.print(io::IO, x::Variable) = show(io, x)
Base.show(io::IO, x::Variable) = print(io, x.name)
function Base.show(io::IO, ::MIME"text/plain", x::Variable)
    known = x.known ? "known" : "unknown"
    print(io, x.name, " (callable ", known, " variable)")
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
function _parse_vars(macroname, known, x)
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
            var_name, expr = _construct_vars(_var.args[1], known, _var.args[2:end])
        else
            var_name, expr = _construct_vars(_var, known, nothing)
        end
        push!(var_names, var_name)
        push!(ex.args, expr)
    end
    push!(ex.args, build_expr(:tuple, var_names))
    return ex
end

function _construct_vars(_var, known, call_args)
    issym  = _var isa Symbol
    isarray = isa(_var, Expr) && _var.head == :ref
    if isarray
        var_name = _var.args[1]
        indices = _var.args[2:end]
        expr = _construct_array_vars(var_name, known, call_args, indices...)
    else
        # Implicit 0-args call
        var_name = _var
        expr = _construct_var(var_name, known, call_args)
    end
    var_name, :($var_name = $expr)
end

function _construct_var(var_name, known, call_args)
    if call_args === nothing
        :(Variable($(Meta.quot(var_name)); known = $known)())
    elseif !isempty(call_args) && call_args[end] == :..
        :(Variable($(Meta.quot(var_name)); known = $known))
    else
        :(Variable($(Meta.quot(var_name)); known = $known)($(call_args...)))
    end
end

function _construct_var(var_name, known, call_args, ind)
    if call_args === nothing
        :(Variable($(Meta.quot(var_name)), $ind...; known = $known)())
    elseif !isempty(call_args) && call_args[end] == :..
        :(Variable($(Meta.quot(var_name)), $ind...; known = $known))
    else
        :(Variable($(Meta.quot(var_name)), $ind...; known = $known)($(call_args...)))
    end
end


function _construct_array_vars(var_name, known, call_args, indices...)
    :(map(Iterators.product($(indices...))) do ind
        $(_construct_var(var_name, known, call_args, :ind))
    end)
end


"""
$(SIGNATURES)

Define one or more unknown variables.
"""
macro variables(xs...)
    esc(_parse_vars(:variables, false, xs))
end

"""
$(SIGNATURES)

Define one or more known variables.
"""
macro parameters(xs...)
    esc(_parse_vars(:parameters, true, xs))
end
