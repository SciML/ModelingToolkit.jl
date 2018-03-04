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

function count_order(x)
    @assert !(x isa Symbol) "The variable $x must have a order of differentiation that is greater or equal to 1!"
    n = 1
    while !(x.args[1] isa Symbol)
        n = n+1
        x = x.args[1]
    end
    n, x.args[1]
end

function _differetial_macro(x)
    ex = Expr(:block)
    lhss = Symbol[]
    for di in x
        @assert di isa Expr && di.args[1] == :~ "@Deriv expects a form that looks like `@Deriv D''~t E'~t`"
        lhs = di.args[2]
        rhs = di.args[3]
        order, lhs = count_order(lhs)
        push!(lhss, lhs)
        expr = :($lhs = Differential($rhs, $order))
        push!(ex.args,  expr)
    end
    push!(ex.args, Expr(:tuple, lhss...))
    ex
end

macro Deriv(x...)
    esc(_differetial_macro(x))
end

export Differential, Derivative, expand_derivatives, @Deriv
