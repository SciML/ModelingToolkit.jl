export Variable, Unknown, Parameter, @Unknown, @Param


struct Variable <: Expression
    name::Symbol
    known::Bool
    dependents::Vector{Variable}
end

Parameter(name, dependents = Variable[]) = Variable(name, true, dependents)
Unknown(name, dependents = Variable[]) = Variable(name, false, dependents)


struct Constant <: Expression
    value::Number
end
Base.get(c::Constant) = c.value


Base.iszero(ex::Expression) = isa(ex, Constant) && iszero(ex.value)
Base.isone(ex::Expression)  = isa(ex, Constant) && isone(ex.value)


# Variables use isequal for equality since == is an Operation
Base.isequal(x::Variable, y::Variable) = (x.name, x.known) == (y.name, y.known)
Base.isequal(::Variable, ::Number) = false
Base.isequal(::Number, ::Variable) = false
Base.isequal(::Variable, ::Constant) = false
Base.isequal(::Constant, ::Variable) = false
Base.isequal(c::Constant, n::Number) = c.value == n
Base.isequal(n::Number, c::Constant) = c.value == n
Base.isequal(a::Constant, b::Constant) = a.value == b.value

function Base.convert(::Type{Expr}, x::Variable)
    x.known               || return x.name
    isempty(x.dependents) && return x.name
    return :($(x.name)($(convert.(Expr, x.dependents)...)))
end
Base.convert(::Type{Expr}, c::Constant) = c.value

Base.show(io::IO, x::Variable) = print(io, x.name)

# Build variables more easily
function _parse_vars(macroname, fun, x)
    ex = Expr(:block)
    var_names = Symbol[]
    # if parsing things in the form of
    # begin
    #     x
    #     y
    #     z
    # end
    x = flatten_expr!(x)
    for _var in x
        iscall = isa(_var, Expr) && _var.head == :call
        issym    = _var isa Symbol
        @assert iscall || issym "@$macroname expects a tuple of expressions (`@$macroname x y z(t)`)"

        if iscall
            dependents = :(Variable[$(_var.args[2:end]...)])
            var_name = _var.args[1]
        else
            dependents = Variable[]
            var_name = _var
        end

        push!(var_names, var_name)
        expr = :($var_name = $fun($(Meta.quot(var_name)), $dependents))
        push!(ex.args, expr)
    end
    push!(ex.args, build_expr(:tuple, var_names))
    return ex
end
macro Unknown(xs...)
    esc(_parse_vars(:Unknown, Unknown, xs))
end
macro Param(xs...)
    esc(_parse_vars(:Param, Parameter, xs))
end
