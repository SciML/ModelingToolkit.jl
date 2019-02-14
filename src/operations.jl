struct Operation <: Expression
    op::Function
    args::Vector{Expression}
end

# Recursive ==
function Base.isequal(x::Operation,y::Operation)
    x.op == y.op && length(x.args) == length(y.args) && all(isequal.(x.args,y.args))
end
Base.isequal(::Operation, ::Number   ) = false
Base.isequal(::Number   , ::Operation) = false
Base.isequal(::Operation, ::Variable ) = false
Base.isequal(::Variable , ::Operation) = false
Base.isequal(::Operation, ::Constant ) = false
Base.isequal(::Constant , ::Operation) = false

Base.convert(::Type{Expr}, O::Operation) =
    build_expr(:call, Any[Symbol(O.op); convert.(Expr, O.args)])
Base.show(io::IO, O::Operation) = print(io, convert(Expr, O))

# For inv
Base.convert(::Type{Operation}, x::Number) = Operation(identity, Expression[Constant(x)])
Base.convert(::Type{Operation}, x::Operation) = x
Base.convert(::Type{Operation}, x::Expression) = Operation(identity, Expression[x])
Operation(x) = convert(Operation, x)
