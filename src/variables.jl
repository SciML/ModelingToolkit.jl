export Variable, @variables, @parameters


struct Variable <: Function
    name::Symbol
    known::Bool
    Variable(name; known = false) = new(name, known)
end
(x::Variable)(args...) = Operation(x, collect(Expression, args))

Base.isequal(x::Variable, y::Variable) = (x.name, x.known) == (y.name, y.known)
Base.show(io::IO, x::Variable) = print(io, x.name)


struct Constant <: Expression
    value::Number
end
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
            expr = :($var_name = $Variable($(Meta.quot(var_name)); known = $known)($(_var.args[2:end]...)))
        else
            var_name = _var
            expr = :($var_name = $Variable($(Meta.quot(var_name)); known = $known))
        end

        push!(var_names, var_name)
        push!(ex.args, expr)
    end
    push!(ex.args, build_expr(:tuple, var_names))
    return ex
end
macro variables(xs...)
    esc(_parse_vars(:variables, false, xs))
end
macro parameters(xs...)
    esc(_parse_vars(:parameters, true, xs))
end
