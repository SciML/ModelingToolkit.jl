using SymbolicUtils: Term, symtype, Sym, simplify

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
    x
    Differential(x) = new(value(x))
end
(D::Differential)(x) = Term{symtype(x)}(D, [x])
(D::Differential)(x::Num) = Num(D(value(x)))

Base.show(io::IO, D::Differential) = print(io, "(D'~", D.x, ")")

Base.:(==)(D1::Differential, D2::Differential) = isequal(D1.x, D2.x)

_isfalse(occ::Bool) = occ === false
_isfalse(occ::Term) = _isfalse(occ.op)

function occursin_info(x, expr::Term)
    if isequal(x, expr)
        true
    else
        args = map(a->occursin_info(x, a), expr.args)
        if all(_isfalse, args)
            return false
        end
        Term{Real}(true, args)
    end
end
function occursin_info(x, expr::Sym)
    isequal(x, expr)
end

hasderiv(O::Term) = O.op isa Differential || any(hasderiv, O.args)
hasderiv(O) = false

occursin_info(x, y) = false
"""
$(SIGNATURES)

TODO
"""
function expand_derivatives(O::Term, simplify=true; occurances=nothing)
    if isa(O.op, Differential)
        @assert length(O.args) == 1
        arg = expand_derivatives(O.args[1], false)

        if occurances == nothing
            occurances = occursin_info(O.op.x, arg)
        end

        _isfalse(occurances) && return 0
        occurances isa Bool && return 1 # means it's a `true`

        (D, o) = (O.op, arg)

        if !isa(o, Term)
            return O # Cannot expand
        elseif isa(o.op, Sym)
            return O # Cannot expand
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
        exprs = []
        c = 0

        for i in 1:l
            t2 = expand_derivatives(D(o.args[i]),false, occurances=occurances.args[i])

            x = if _iszero(t2)
                t2
            elseif _isone(t2)
                d = derivative_idx(o, i)
                d isa NoDeriv ? D(o) : d
            else
                t1 = derivative_idx(o, i)
                t1 = t1 isa NoDeriv ? D(o) : t1
                make_operation(*, [t1, t2])
            end

            if _iszero(x)
                continue
            elseif x isa Symbolic
                push!(exprs, x)
            else
                c += x
            end
        end

        if isempty(exprs)
            return c
        elseif length(exprs) == 1
            term = (simplify ? SymbolicUtils.simplify(exprs[1]) : exprs[1])
            return _iszero(c) ? term : c + term
        else
            x = make_operation(+, !_iszero(c) ? vcat(c, exprs) : exprs)
            return simplify ? SymbolicUtils.simplify(x) : x
        end
    elseif !hasderiv(O)
        return O
    else
        args = map(a->expand_derivatives(a, false), O.args)
        O1 = make_operation(O.op, args)
        return simplify ? SymbolicUtils.simplify(O1) : O1
    end
end

function expand_derivatives(n::Num, simplify=true; occurances=nothing)
    Num(expand_derivatives(value(n), simplify; occurances=occurances))
end

_iszero(x) = false
_isone(x) = false

expand_derivatives(x, simplify=true;occurances=nothing) = x

# Don't specialize on the function here
"""
$(SIGNATURES)

Calculate the derivative of the op `O` with respect to its argument with index
`idx`.

# Examples

```jldoctest label1
julia> using ModelingToolkit

julia> @variables x y;

julia> ModelingToolkit.derivative_idx(sin(x), 1)
cos(x())
```

Note that the function does not recurse into the operation's arguments, i.e., the
chain rule is not applied:

```jldoctest label1
julia> myop = sin(x) * y^2
sin(x()) * y() ^ 2

julia> typeof(myop.op)  # Op is multiplication function
typeof(*)

julia> ModelingToolkit.derivative_idx(myop, 1)  # wrt. sin(x)
y() ^ 2

julia> ModelingToolkit.derivative_idx(myop, 2)  # wrt. y^2
sin(x())
```
"""
derivative_idx(O::Any, ::Any) = 0
derivative_idx(O::Term, idx) = derivative(O.op, (O.args...,), Val(idx))

# Indicate that no derivative is defined.
struct NoDeriv
end
derivative(f, args, v) = NoDeriv()

# Pre-defined derivatives
import DiffRules
for (modu, fun, arity) ∈ DiffRules.diffrules()
    fun in [:*, :+, :abs, :mod, :rem, :max, :min] && continue # special
    for i ∈ 1:arity

        expr = if arity == 1
            DiffRules.diffrule(modu, fun, :(args[1]))
        else
            DiffRules.diffrule(modu, fun, ntuple(k->:(args[$k]), arity)...)[i]
        end
        @eval derivative(::typeof($modu.$fun), args::NTuple{$arity,Any}, ::Val{$i}) = $expr
    end
end

derivative(::typeof(+), args::NTuple{N,Any}, ::Val) where {N} = 1
derivative(::typeof(*), args::NTuple{N,Any}, ::Val{i}) where {N,i} = make_operation(*, deleteat!(collect(args), i))
derivative(::typeof(one), args::Tuple{<:Any}, ::Val) = 0

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
        expr = :($lhs = $_repeat_apply(Differential($value($rhs)), $order))
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
