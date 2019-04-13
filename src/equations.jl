export Equation


struct Equation
    lhs::Expression
    rhs::Expression
end
Base.:(==)(a::Equation, b::Equation) = isequal((a.lhs, a.rhs), (b.lhs, b.rhs))

Base.:~(lhs::Expression, rhs::Expression) = Equation(lhs, rhs)
Base.:~(lhs::Expression, rhs::Number    ) = Equation(lhs, rhs)
Base.:~(lhs::Number    , rhs::Expression) = Equation(lhs, rhs)

function find_parameter_calls(O::Operation,p_calls=Variable[])
  if O.op isa Variable && O.op.known && !isempty(O.args)
      push!(p_calls,O.op)
      find_parameter_calls.(O.args,(p_calls,))
  else
      find_parameter_calls.(O.args,(p_calls,))
  end
  p_calls
end
find_parameter_calls(O,p_calls) = nothing
