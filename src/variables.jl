mutable struct Variable <: Expression
    name::Symbol
    subtype::Symbol
    diff::Union{Function,Nothing}  # FIXME
    dependents::Vector{Variable}
end

Variable(name; subtype::Symbol=:Variable, dependents::Vector{Variable} = Variable[]) =
    Variable(name, subtype, nothing, dependents)
Variable(name, args...; kwargs...) = Variable(name, args...; subtype=:Variable, kwargs...)

Parameter(name,args...;kwargs...) = Variable(name,args...;subtype=:Parameter,kwargs...)
IndependentVariable(name,args...;kwargs...) = Variable(name,args...;subtype=:IndependentVariable,kwargs...)

function DependentVariable(name,args...;dependents = [],kwargs...)
    @assert !isempty(dependents)
    Variable(name,args...;subtype=:DependentVariable,dependents=dependents,kwargs...)
end

export Variable,Parameter,Constant,DependentVariable,IndependentVariable,
       @Var, @Param, @Const, @DVar, @IVar


Base.copy(x::Variable) = Variable(x.name, x.subtype, x.diff, x.dependents)


struct Constant <: Expression
    value::Number
end
Base.get(c::Constant) = c.value


Base.iszero(ex::Expression) = isa(ex, Constant) && iszero(ex.value)
Base.isone(ex::Expression)  = isa(ex, Constant) && isone(ex.value)


# Variables use isequal for equality since == is an Operation
Base.:(==)(x::Variable, y::Variable) = (x.name, x.subtype, x.diff) == (y.name, y.subtype, y.diff)
Base.:(==)(::Variable, ::Number) = false
Base.:(==)(::Number, ::Variable) = false
Base.:(==)(::Variable, ::Constant) = false
Base.:(==)(::Constant, ::Variable) = false
Base.:(==)(c::Constant, n::Number) = c.value == n
Base.:(==)(n::Number, c::Constant) = c.value == n
Base.:(==)(a::Constant, b::Constant) = a.value == b.value

function Base.convert(::Type{Expr}, x::Variable)
    x.diff === nothing && return x.name
    return Symbol("$(x.name)_$(x.diff.x.name)")
end
Base.convert(::Type{Expr}, c::Constant) = c.value

function Base.show(io::IO, x::Variable)
    print(io, x.subtype, '(', x.name, ')')
    x.diff === nothing || print(io, ", diff = ", x.diff)
end

# Build variables more easily
function _parse_vars(macroname, fun, x)
    ex = Expr(:block)
    lhss = Symbol[]
    # if parsing things in the form of
    # begin
    #     x
    #     y
    #     z
    # end
    x = flatten_expr!(x)
    for _var in x
        iscall = typeof(_var) <: Expr && _var.head == :call
        issym    = _var isa Symbol
        @assert iscall || issym "@$macroname expects a tuple of expressions!\nE.g. `@$macroname x y z`"

        if iscall
            dependents = :([$(_var.args[2:end]...)])
            lhs = _var.args[1]
        else
            dependents = Variable[]
            lhs = _var
        end

        push!(lhss, lhs)
        expr = :( $lhs = $fun( Symbol($(String(lhs))) ,
                  dependents = $dependents))
        push!(ex.args, expr)
    end
    push!(ex.args, build_expr(:tuple, lhss))
    return ex
end

for funs in ((:DVar, :DependentVariable), (:IVar, :IndependentVariable),
             (:Var, :Variable),
             (:Param, :Parameter))
    @eval begin
        macro ($(funs[1]))(x...)
            esc(_parse_vars(String($funs[1]), $funs[2], x))
        end
    end
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
