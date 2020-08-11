"""
$(TYPEDEF)

Represents a differential operator.

# Fields
$(FIELDS)

# Examples

```jldoctest
julia> using ModelingToolkit

julia> @variables x y;

julia> D = Differential(x)
(D'~x())

julia> D(y)  # Differentiate y wrt. x
(D'~x())(y())
```
"""
struct Differential <: Function
    """The variable or expression to differentiate with respect to."""
    x::Expression
end
(D::Differential)(x) = Operation(D, Expression[x])

Base.show(io::IO, D::Differential) = print(io, "(D'~", D.x, ")")
Base.convert(::Type{Expr}, D::Differential) = D

Base.:(==)(D1::Differential, D2::Differential) = isequal(D1.x, D2.x)

_isfalse(occ::Constant) = occ.value === false
_isfalse(occ::Operation) = _isfalse(occ.op)

function occursin_info(x, expr::Operation)
    if isequal(x, expr)
        Constant(true)
    else
        args = map(a->occursin_info(x, a), expr.args)
        if all(_isfalse, args)
            return Constant(false)
        end
        Operation(Constant(true), args)
    end
end

hasderiv(O::Operation) = O.op isa Differential || any(hasderiv, O.args)
hasderiv(O) = false

occursin_info(x, y) = Constant(false)
"""
$(SIGNATURES)

TODO
"""
function expand_derivatives(O::Operation, simplify=true; occurances=nothing)
    if isa(O.op, Differential)
        @assert length(O.args) == 1
        arg = expand_derivatives(O.args[1], false)

        if occurances == nothing
            occurances = occursin_info(O.op.x, arg)
        end

        _isfalse(occurances) && return Constant(0)
        occurances isa Constant && return Constant(1) # means it's a Constant(true)

        (D, o) = (O.op, arg)

        if o isa Constant
            return Constant(0)
        elseif isequal(o, D.x)
            return Constant(1)
        elseif !isa(o, Operation)
            return O
        elseif isa(o.op, Variable)
            return O # This means D.x occurs in o, but o is symbolic function call
        elseif isa(o.op, Differential)
            # The recursive expand_derivatives was not able to remove
            # a nested Differential. We can attempt to differentiate the
            # inner expression wrt to the outer iv. And leave the
            # unexpandable Differential outside.
            if isequal(o.op.x, D.x)
                return O
            else
                return expand_derivatives(o.op(expand_derivatives(D(o.args[1]), false)), simplify)
            end
        end

        l = length(o.args)
        exprs = Expression[]
        c = 0

        for i in 1:l
            t2 = expand_derivatives(D(o.args[i]),false, occurances=occurances.args[i])

            x = if _iszero(t2)
                t2
            elseif _isone(t2)
                derivative(o, i)
            else
                t1 = derivative(o, i)
                make_operation(*, Expression[t1, t2])
            end

            if _iszero(x)
                continue
            elseif x isa Expression
                push!(exprs, x)
            elseif x isa Constant
                c += x.value
            else
                c += x
            end
        end

        if isempty(exprs)
            return Constant(c)
        elseif length(exprs) == 1
            return simplify ? ModelingToolkit.simplify(exprs[1]) : exprs[1]
        else
            x = make_operation(+, !iszero(c) ? vcat(c, exprs) : exprs)
            return simplify ? ModelingToolkit.simplify(x) : x
        end
    elseif !hasderiv(O)
        return O
    else
        args = map(a->expand_derivatives(a, false), O.args)
        O1 = make_operation(O.op, args)
        return simplify ? ModelingToolkit.simplify(O1) : O1
    end
end
_iszero(x::Constant) = iszero(x.value)
_isone(x::Constant) = isone(x.value)
_iszero(x) = false
_isone(x) = false

expand_derivatives(x,args...;occurances=nothing) = x

# Don't specialize on the function here
"""
$(SIGNATURES)

Calculate the derivative of the op `O` with respect to its argument with index
`idx`.

# Examples

```jldoctest label1
julia> using ModelingToolkit

julia> @variables x y;

julia> ModelingToolkit.derivative(sin(x), 1)
cos(x())
```

Note that the function does not recurse into the operation's arguments, i.e., the
chain rule is not applied:

```jldoctest label1
julia> myop = sin(x) * y^2
sin(x()) * y() ^ 2

julia> typeof(myop.op)  # Op is multiplication function
typeof(*)

julia> ModelingToolkit.derivative(myop, 1)  # wrt. sin(x)
y() ^ 2

julia> ModelingToolkit.derivative(myop, 2)  # wrt. y^2
sin(x())
```
"""
derivative(O::Operation, idx) = derivative(O.op, (O.args...,), Val(idx))
derivative(O::Constant, ::Any) = Constant(0)

# Pre-defined derivatives
import DiffRules
for (modu, fun, arity) ∈ DiffRules.diffrules()
    fun in [:*, :+] && continue # special
    for i ∈ 1:arity

        expr = if arity == 1
            DiffRules.diffrule(modu, fun, :(args[1]))
        else
            DiffRules.diffrule(modu, fun, ntuple(k->:(args[$k]), arity)...)[i]
        end
        @eval derivative(::typeof($modu.$fun), args::NTuple{$arity,Any}, ::Val{$i}) = (x = $expr; !(x isa Expression || x isa Constant) ? Constant(x) : x)
    end
end

derivative(::typeof(+), args::NTuple{N,Any}, ::Val) where {N} = Constant(1)
derivative(::typeof(*), args::NTuple{N,Any}, ::Val{i}) where {N,i} = make_operation(*, deleteat!(collect(args), i))
derivative(::typeof(one), args::Tuple{<:Any}, ::Val) = Constant(0)

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

"""
$(SIGNATURES)

Define one or more differentials.

# Examples

```jldoctest
julia> using ModelingToolkit

julia> @variables x y z;

julia> @derivatives Dx'~x Dy'~y  # Create differentials wrt. x and y
((D'~x()), (D'~y()))

julia> Dx(z)  # Differentiate z wrt. x
(D'~x())(z())

julia> Dy(z)  # Differentiate z wrt. y
(D'~y())(z())
```
"""
macro derivatives(x...)
    esc(_differential_macro(x))
end
