export Equation


mutable struct Equation
    lhs::Expression
    rhs::Expression
end
Base.broadcastable(eq::Equation) = Ref(eq)

Base.:~(lhs::Expression, rhs::Expression) = Equation(lhs, rhs)
Base.:~(lhs::Expression, rhs::Number    ) = Equation(lhs, rhs)
Base.:~(lhs::Number    , rhs::Expression) = Equation(lhs, rhs)
