# Parameterize by T so that way it can be Vector{Expression} which is defined after
struct Operation{T}
    op::Function
    args::Vector{T}
end

# Operations use isequal for equality since == is an Operation
function Base.isequal(x::Operation,y::Operation)
    x.op == y.op && all(isequal.(x.args,y.args))
end

function Base.Expr(O::Operation)
    :($(Symbol(O.op))($((Expr(x) for x in O.args)...)))
end

Base.show(io::IO,O::Operation) = print(io,string(Expr(O)))
