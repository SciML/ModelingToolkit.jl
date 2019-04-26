export Variable, @variables, @parameters


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
        issym    = _var isa Symbol
        @assert iscall || issym "@$macroname expects a tuple of expressions or an expression of a tuple (`@$macroname x y z(t)` or `@$macroname x, y, z(t)`)"

        if iscall
            var_name = _var.args[1]
            if _var.args[end] == :..
                expr = :($var_name = $Variable($(Meta.quot(var_name)); known = $known))
            else
                expr = :($var_name = $Variable($(Meta.quot(var_name)); known = $known)($(_var.args[2:end]...)))
            end
        else
            # Implicit 0-args call
            var_name = _var
            expr = :($var_name = $Variable($(Meta.quot(var_name)); known = $known)())
        end

        push!(var_names, var_name)
        push!(ex.args, expr)
    end
    push!(ex.args, build_expr(:tuple, var_names))
    return ex
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
