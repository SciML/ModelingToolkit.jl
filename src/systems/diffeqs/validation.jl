isintermediate(eq::Equation) = !(isa(eq.lhs, Operation) && isa(eq.lhs.op, Differential))

struct ODEExpr  # dⁿx/dtⁿ = rhs
    x::Variable
    n::Int
    rhs::Expression
end
function Base.convert(::Type{ODEExpr},eq::Equation)
    isintermediate(eq) && throw(ArgumentError("intermediate equation received"))
    (x, t, n) = flatten_differential(eq.lhs)
    (isa(t, Operation) && isa(t.op, Variable) && isempty(t.args)) ||
        throw(ArgumentError("invalid independent variable $t"))
    (isa(x, Operation) && isa(x.op, Variable) && length(x.args) == 1 && isequal(first(x.args), t)) ||
        throw(ArgumentError("invalid dependent variable $x"))
    return t.op, ODEExpr(x.op, n, eq.rhs)
end
Base.:(==)(a::ODEExpr, b::ODEExpr) = isequal((a.x, a.n, a.rhs), (b.x, b.n, b.rhs))
