export Differential, expand_derivatives, @derivatives


struct Differential <: Function
    x::Expression
end

Base.show(io::IO, D::Differential) = print(io, "(D'~", D.x, ")")
Base.convert(::Type{Expr}, D::Differential) = D

(D::Differential)(x::Operation) = Operation(D, Expression[x])
function (D::Differential)(x::Variable)
    D.x === x             && return Constant(1)
    has_dependent(x, D.x) || return Constant(0)
    return Operation(D, Expression[x])
end
(::Differential)(::Any) = Constant(0)
Base.:(==)(D1::Differential, D2::Differential) = isequal(D1.x, D2.x)

function expand_derivatives(O::Operation)
    @. O.args = expand_derivatives(O.args)

    if O.op isa Differential
        D = O.op
        o = O.args[1]
        isa(o, Operation) || return O
        return simplify_constants(sum(i->derivative(o,i)*expand_derivatives(D(o.args[i])),1:length(o.args)))
    end

    return O
end
expand_derivatives(x) = x

# Don't specialize on the function here
derivative(O::Operation, idx) = derivative(O.op, (O.args...,), Val(idx))

# Pre-defined derivatives
import DiffRules, SpecialFunctions, NaNMath
for (modu, fun, arity) ∈ DiffRules.diffrules()
    for i ∈ 1:arity
        @eval function derivative(::typeof($modu.$fun), args::NTuple{$arity,Any}, ::Val{$i})
            M, f = $(modu, fun)
            partials = DiffRules.diffrule(M, f, args...)
            dx = @static $arity == 1 ? partials : partials[$i]
            convert(Expression, dx)
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

_repeat_apply(f, n) = n == 1 ? f : f ∘ _repeat_apply(f, n-1)
function _differential_macro(x)
    ex = Expr(:block)
    lhss = Symbol[]
    x = x isa Tuple && first(x).head == :tuple ? first(x).args : x # tuple handling
    x = flatten_expr!(x)
    for di in x
        @assert di isa Expr && di.args[1] == :~ "@derivatives expects a form that looks like `@derivatives D''~t E'~t` or `@derivatives (D''~t), (E'~t)`"
        lhs = di.args[2]
        rhs = di.args[3]
        order, lhs = count_order(lhs)
        push!(lhss, lhs)
        expr = :($lhs = $_repeat_apply(Differential($rhs), $order))
        push!(ex.args,  expr)
    end
    push!(ex.args, Expr(:tuple, lhss...))
    ex
end

macro derivatives(x...)
    esc(_differential_macro(x))
end

function calculate_jacobian(eqs,vars)
    Expression[Differential(vars[j])(eqs[i]) for i in 1:length(eqs), j in 1:length(vars)]
end
