# Parameterize by T so that way it can be Vector{Expression} which is defined after
struct Operation <: AbstractOperation
    op::Function
    args::Vector{Expression}
end

# Recursive ==
function Base.:(==)(x::Operation,y::Operation)
    x.op == y.op && all(isequal.(x.args,y.args))
end
Base.:(==)(x::Operation,y::Number) = false
Base.:(==)(x::Number,y::Operation) = false
Base.:(==)(x::Operation,y::Void) = false
Base.:(==)(x::Void,y::Operation) = false
Base.:(==)(x::Variable,y::Operation) = y == Operation(identity,Expression[x])
Base.:(==)(x::Operation,y::Variable) = x == Operation(identity,Expression[y])

Base.convert(::Type{Operation},x::Variable) = Operation(identity,x)


# Don't recurse inversion for Jacobians
Base.inv(x::Operation) = x

function Base.Expr(O::Operation)
    Expr(:call, Symbol(O.op), Expr.(O.args)...)
end

Base.show(io::IO,O::Operation) = print(io,string(Expr(O)))

# Bigger printing
# Is there a way to just not have this as the default?
function Base.parse(::Type{Operation},ex::Expr)
    f = ex.args[1]
    operands = ex.args[2:end]
    if ex.head == :call && any(x -> x isa Expr, ex.args)
        args = Expression[parse(Operation,o) for o in operands]
        parse(Operation, f, args)
    else
        parse(Operation, f, Expression[parse(Operation,o) for o in operands])
    end
end
Base.parse(::Type{Operation},x::Expression) = x
Base.parse(::Type{Operation},sym::Symbol,args) = Operation(eval(sym), args)
Base.parse(::Type{Operation},x::Union{Symbol, Number}) = x

"""
find_replace(O::Operation,x::Variable,y::Expression)

Finds the variable `x` in Operation `O` and replaces it with the Expression `y`
"""
function find_replace!(O::Operation,x::Variable,y::Expression)
    for i in eachindex(O.args)
        if isequal(O.args[i],x)
            O.args[i] = y
        elseif typeof(O.args[i]) <: Operation
            find_replace!(O.args[i],x,y)
        end
    end
end
