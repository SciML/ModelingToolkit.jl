"""
$(TYPEDEF)

An equality relationship between two expressions.

# Fields
$(FIELDS)
"""
struct Equation
    """The expression on the left-hand side of the equation."""
    lhs
    """The expression on the right-hand side of the equation."""
    rhs
end
Base.:(==)(a::Equation, b::Equation) = all(isequal.((a.lhs, a.rhs), (b.lhs, b.rhs)))
Base.hash(a::Equation, salt::UInt) = hash(a.lhs, hash(a.rhs, salt))

"""
$(TYPEDSIGNATURES)

Create an [`Equation`](@ref) out of two [`Expression`](@ref) instances, or an
`Expression` and a `Number`.

# Examples

```jldoctest
julia> using ModelingToolkit

julia> @variables x y;

julia> x ~ y
Equation(x(), y())

julia> x - y ~ 0
Equation(x() - y(), 0)
```
"""
Base.:~(lhs::Num, rhs::Num) = Equation(value(lhs), value(rhs))
Base.:~(lhs::Num, rhs::Number    ) = Equation(value(lhs), value(rhs))
Base.:~(lhs::Number    , rhs::Num) = Equation(value(lhs), value(rhs))
Base.:~(lhs::Symbolic, rhs::Symbolic) = Equation(value(lhs), value(rhs))
Base.:~(lhs::Symbolic, rhs::Any    ) = Equation(value(lhs), value(rhs))
Base.:~(lhs::Any, rhs::Symbolic    ) = Equation(value(lhs), value(rhs))
Base.:~(lhs::Number    , rhs::Num) = Equation(value(lhs), value(rhs))

struct ConstrainedEquation
  constraints
  eq
end


function _eq_unordered(a, b)
    length(a) === length(b) || return false
    n = length(a)
    idxs = Set(1:n)
    for x ∈ a
        idx = findfirst(isequal(x), b)
        idx === nothing && return false
        idx ∈ idxs      || return false
        delete!(idxs, idx)
    end
    return true
end

function expand_derivatives(eq::Equation, simplify=true)
    return Equation(expand_derivatives(eq.lhs, simplify), expand_derivatives(eq.rhs, simplify))
end
