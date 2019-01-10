struct Differential <: Function
    x::Expression
    order::Int
end
Differential(x) = Differential(x,1)

Base.show(io::IO, D::Differential) = print(io,"($(D.x),$(D.order))")
Base.convert(::Type{Expr}, D::Differential) = D

function Derivative end
(D::Differential)(x) = @term(D(x))
function (D::Differential)(x::Variable)
    D.x === x             && return @term(1)
    has_dependent(x, D.x) || return @term(0)
    return Term(Variable(x, D))
end
Base.:(==)(D1::Differential, D2::Differential) = D1.order == D2.order && D1.x == D2.x

Variable(x::Variable, D::Differential) = Variable(x.name,x.value,x.value_type,
                x.subtype,D,x.dependents,x.description,x.flow,x.domain,
                x.size,x.context)

function expand_derivatives(t::Term)
    is_branch(t) || return t

    head = root(t)
    args = expand_derivatives.(children(t))

    if head === :call && length(args) === 2 && isa(root(args[1]), Differential)
        D = root(args[1])::Differential
        o = args[2]
        o_root = root(o)
        # @info "test" o_root typeof(o_root)
        isa(o_root, Variable) && return D(o_root)
        o_args = children(o)[2:end]
        # @info "expand_derivatives" D o o_root o_args
        args = map(eachindex(o_args)) do i
            @term($(Derivative(o, i)) * $(expand_derivatives(D(o_args[i])))).x
        end

        if length(args) == 1
            ex = first(args)
        else
            ex = Expr(:call)
            push!(ex.args, +)
            append!(ex.args, args)
        end

        return convert(Term, ex) |> simplify_constants
    end

    return t
end
expand_derivatives(x::Variable) = x

# Don't specialize on the function here
function Derivative(t::Term, idx)
    # This calls the Derivative dispatch from the user or pre-defined code
    fn, args = unpack(t)
    Derivative(fn, (args...,), Val(idx))
end
Derivative(fn, args, idx) = Derivative(fn, (args...,), idx)

# Pre-defined derivatives
import DiffRules, SpecialFunctions, NaNMath
for (modu, fun, arity) ∈ DiffRules.diffrules()
    for i ∈ 1:arity
        mc = :(@eval @term)
        push!(mc.args[3].args, Expr(:$, :dx))
        @eval function Derivative(::typeof($modu.$fun), args::NTuple{$arity,Any}, ::Val{$i})
            M, f = $(modu, fun)
            partials = DiffRules.diffrule(M, f, args...)
            dx = @static $arity == 1 ? partials : partials[$i]
            $mc
        end
    end
end
# _eval(t)

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
    Term[Differential(vars[j])(eqs[i]) for i ∈ eachindex(eqs), j ∈ eachindex(vars)]
end

export Differential, Derivative, expand_derivatives, @Deriv, calculate_jacobian
