struct Variable <: Expression
    name::Symbol
    subtype::Symbol
    dependents::Vector{Variable}
end

Parameter(name; dependents = Variable[]) = Variable(name, :Parameter, dependents)
Unknown(name; dependents = Variable[]) = Variable(name, :Unknown, dependents)

export Variable, Unknown, Parameter, Constant, @Unknown, @Param, @Const


struct Constant <: Expression
    value::Number
end
Base.get(c::Constant) = c.value


Base.iszero(ex::Expression) = isa(ex, Constant) && iszero(ex.value)
Base.isone(ex::Expression)  = isa(ex, Constant) && isone(ex.value)


# Variables use isequal for equality since == is an Operation
Base.:(==)(x::Variable, y::Variable) = (x.name, x.subtype) == (y.name, y.subtype)
Base.:(==)(::Variable, ::Number) = false
Base.:(==)(::Number, ::Variable) = false
Base.:(==)(::Variable, ::Constant) = false
Base.:(==)(::Constant, ::Variable) = false
Base.:(==)(c::Constant, n::Number) = c.value == n
Base.:(==)(n::Number, c::Constant) = c.value == n
Base.:(==)(a::Constant, b::Constant) = a.value == b.value

Base.convert(::Type{Expr}, x::Variable) = x.name
Base.convert(::Type{Expr}, c::Constant) = c.value

Base.show(io::IO, x::Variable) = print(io, x.subtype, '(', x.name, ')')

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
        expr = :($var_name = $fun($(Meta.quot(var_name)), dependents = $dependents))
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

function _const_assign(x)
    ex = Expr(:block)
    lhss = Symbol[]
    x = flatten_expr!(x)
    for eq in x
        @assert eq isa Expr && eq.head == :(=) "@Const expects a tuple of assignments!\nE.g. `@Const D=t W=g`"
        lhs = eq.args[1]
        push!(lhss, lhs)
        rhs = eq.args[2]
        expr = :($lhs = Constant($rhs))
        push!(ex.args,  expr)
    end
    push!(ex.args, Expr(:tuple, lhss...))
    ex
end

macro Const(x...)
    esc(_const_assign(x))
end
