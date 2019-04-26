export Differential, expand_derivatives, @derivatives


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

"""
$(SIGNATURES)

TODO
"""
function expand_derivatives(O::Operation)
    @. O.args = expand_derivatives(O.args)

    if isa(O.op, Differential)
        (D, o) = (O.op, O.args[1])

        isequal(o, D.x)     && return Constant(1)
        occursin(D.x, o)    || return Constant(0)
        isa(o, Operation)   || return O
        isa(o.op, Variable) && return O

        return sum(1:length(o.args)) do i
            derivative(o, i) * expand_derivatives(D(o.args[i]))
        end |> simplify_constants
    end

    return O
end
expand_derivatives(x) = x

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

Note that the function does not recurse into the operation's arguments, i.e. the
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


function calculate_jacobian(eqs, dvs)
    Expression[Differential(dv)(eq) for eq ∈ eqs, dv ∈ dvs]
end
