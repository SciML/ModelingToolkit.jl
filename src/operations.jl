struct Operation <: Expression
    op::Function
    args::Vector{Expression}
end

# Recursive ==
function Base.:(==)(x::Operation,y::Operation)
    x.op == y.op && length(x.args) == length(y.args) && all(isequal.(x.args,y.args))
end
Base.:(==)(::Operation, ::Number   ) = false
Base.:(==)(::Number   , ::Operation) = false
Base.:(==)(::Operation, ::Variable ) = false
Base.:(==)(::Variable , ::Operation) = false
Base.:(==)(::Operation, ::Constant ) = false
Base.:(==)(::Constant , ::Operation) = false

Base.convert(::Type{Expr}, O::Operation) =
    build_expr(:call, Any[Symbol(O.op); convert.(Expr, O.args)])
Base.show(io::IO, O::Operation) = print(io, convert(Expr, O))

# For inv
Base.convert(::Type{Operation}, x::Number) = Operation(identity, Expression[Constant(x)])
Base.convert(::Type{Operation}, x::Operation) = x
Base.convert(::Type{Operation}, x::Expression) = Operation(identity, Expression[x])
Operation(x) = convert(Operation, x)
