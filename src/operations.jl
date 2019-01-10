# Parameterize by T so that way it can be Vector{Expression} which is defined after
struct Operation <: AbstractOperation
    op::Function
    args::Vector{Expression}
end

# Recursive ==
function Base.:(==)(x::Operation,y::Operation)
    x.op == y.op && length(x.args) == length(y.args) && all(isequal.(x.args,y.args))
end
Base.:(==)(x::Operation,y::Number) = false
Base.:(==)(x::Number,y::Operation) = false
Base.:(==)(x::Operation,y::Nothing) = false
Base.:(==)(x::Nothing,y::Operation) = false
Base.:(==)(x::Variable,y::Operation) = false
Base.:(==)(x::Operation,y::Variable) = false

Base.convert(::Type{Expr}, O::Operation) = Expr(:call, Symbol(O.op), convert.(Expr, O.args)...)
Base.show(io::IO, O::Operation) = print(io, convert(Expr, O))


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

# For inv
Base.convert(::Type{Operation},x::Int) = Operation(identity,Expression[Constant(x)])
Base.convert(::Type{Operation},x::Bool) = Operation(identity,Expression[Constant(x)])
Base.convert(::Type{Operation},x::Variable) = Operation(identity,Expression[x])
Operation(x) = convert(Operation,x)
Operation(x::Operation) = x
