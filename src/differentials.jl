struct Differential <: Expression
    x::Expression
    order::Int
end
Differential(x) = Differential(x,1)

Base.show(io::IO, D::Differential) = print(io,"($(D.x),$(D.order))")
Base.Expr(D::Differential) = :($(Symbol("D_$(D.x.name)_$(D.order)")))

function Derivative end
(D::Differential)(x::Operation) = Operation(Derivative,Expression[x,D])
function (D::Differential)(x::Variable)
    D.x === x             && return Constant(1)
    has_dependent(x, D.x) || return Constant(0)
    return Variable(x,D)
end
Base.:(==)(D1::Differential, D2::Differential) = D1.order == D2.order && D1.x == D2.x

Variable(x::Variable,D::Differential) = Variable(x.name,x.value,x.value_type,
                x.subtype,D,x.dependents,x.description,x.flow,x.domain,
                x.size,x.context)

function expand_derivatives(O::Operation)
    @. O.args = expand_derivatives(O.args)

    if O.op == Derivative
        D = O.args[2]
        o = O.args[1]
        return simplify_constants(sum(i->Derivative(o,i)*expand_derivatives(D(o.args[i])),1:length(o.args)))
    end

    return O
end
expand_derivatives(x::Variable) = x
expand_derivatives(D::Differential) = D

# Don't specialize on the function here
function Derivative(O::Operation,idx)
    # This calls the Derivative dispatch from the user or pre-defined code
    Derivative(O.op, O.args, Val(idx))
end
Derivative(op, args, idx) = Derivative(op, (args...,), idx)

# Pre-defined derivatives
import DiffRules, SpecialFunctions, NaNMath
for (modu, fun, arity) ∈ DiffRules.diffrules()
    for i ∈ 1:arity
        @eval function Derivative(::typeof($modu.$fun), args::NTuple{$arity,Any}, ::Val{$i})
            M, f = $(modu, fun)
            partials = DiffRules.diffrule(M, f, args...)
            dx = @static $arity == 1 ? partials : partials[$i]
            parse(Operation,dx)
        end
    end
end

function count_order(x)
    @assert !(x isa Symbol) "The variable $x must have an order of differentiation that is greater or equal to 1!"
    n = 1
    while !(x.args[1] isa Symbol)
        n = n+1
        x = x.args[1]
    end
    n, x.args[1]
end

function _differential_macro(x)
    ex = Expr(:block)
    lhss = Symbol[]
    x = flatten_expr!(x)
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
    esc(_differential_macro(x))
end

function calculate_jacobian(eqs,vars)
    Expression[Differential(vars[j])(eqs[i]) for i in 1:length(eqs), j in 1:length(vars)]
end

export Differential, Derivative, expand_derivatives, @Deriv, calculate_jacobian
