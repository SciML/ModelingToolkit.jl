export Equation


"""
$(TYPEDEF)

An equality relationship between two expressions.

# Fields
$(FIELDS)
"""
struct Equation
    """The expression on the left hand side of the equation."""
    lhs::Expression
    """The expression on the right hand side of the equation."""
    rhs::Expression
end
Base.:(==)(a::Equation, b::Equation) = isequal((a.lhs, a.rhs), (b.lhs, b.rhs))

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
Equation(x() - y(), ModelingToolkit.Constant(0))
```
"""
Base.:~(lhs::Expression, rhs::Expression) = Equation(lhs, rhs)
Base.:~(lhs::Expression, rhs::Number    ) = Equation(lhs, rhs)
Base.:~(lhs::Number    , rhs::Expression) = Equation(lhs, rhs)
