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
    Expr(:call, Symbol(O.op), O.args...)
end

Base.show(io::IO,O::Operation) = print(io,string(Expr(O)))

Operation(sym::Symbol, args) = Operation(eval(sym), args)
Operation(x::Union{Symbol, Number}) = x

function Operation(ex::Expr)
    f = ex.args[1]
    operands = ex.args[2:end]
    if ex.head == :call && any(x -> x isa Expr, ex.args)
        args = [map(Operation, operands)]
        Operation(f, args...)
    else
        Operation(f, operands)
    end
end
