struct Differential <: AbstractOperator
    x::Union{Variable,Operation}
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
        #=
        diff_idxs = find(x->isequal(x,by.x),O.args)
        (diff_idxs != nothing || length(diff_idxs) > 1) && error("Derivatives of multi-argument functions require matching a unique argument.")
        idx = first(diff_idxs)
        =#
        i = 1
        if typeof(O.args[1].args[i]) == typeof(O.args[2].x) && isequal(O.args[1].args[i],O.args[2].x)
            Derivative(O.args[1],i)
        else
            D = Differential(O.args[2].x)
            cr_exp = D*O.args[1].args[i]
            Derivative(O.args[1],i) * expand_derivatives(cr_exp)
        end
    else
        for i in 1:length(O.args)
            O.args[i] = expand_derivatives(O.args[i])
        end
    end
end

# Don't specialize on the function here
function Derivative(O::Operation,idx)
    Derivative(O.op,O.args,Val{idx})
end

## Pre-defined derivatives
function Derivative(::typeof(sin),args,::Type{Val{1}})
    Operation(Base.cos,args)
end

export Differential, Derivative, expand_derivatives
