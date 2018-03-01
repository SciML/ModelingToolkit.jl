struct Differential <: AbstractOperator
    x::Variable
    order::Int
end
Differential(x) = Differential(x,1)

Base.show(io::IO, D::Differential) = print(io,"($(D.x),$(D.order))")
Base.Expr(D::Differential) = :($(Symbol("D_$(D.x.name)_$(D.order)")))

function Derivative end
Base.:*(D::Differential,x::Operation) = Operation(Derivative,Expression[x,D])
function Base.:*(D::Differential,x::Variable)
    if D.x === x
        Constant(1)
    else
        Variable(x.name,x.subtype,x.value,x.value_type,D)
    end
end

function expand_derivatives(O::Operation)
    if O.op == Derivative
        @assert length(O.args) == 2
        Derivative(O.args[1],O.args[2])
    else
        for i in 1:length(O.args)
            O.args[i] = expand_derivatives(O.args[i])
        end
    end
end

# Don't specialize on the function here
function Derivative(O::Operation,by::Differential)
    diff_idxs = find(x->isequal(x,by.x),O.args)
    @assert diff_idxs != nothing && length(diff_idxs) == 1
    idx = first(diff_idxs)
    Derivative(O.op,O.args,idx)
end

## Pre-defined derivatives
function Derivative(::typeof(sin),args,idx)
    Operation(Base.cos,args)
end

export Differential, Derivative, expand_derivatives
