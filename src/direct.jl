"""
```julia
gradient(O::Expression, vars::AbstractVector{<:Expression}; simplify = true)
```

A helper function for computing the gradient of an expression with respect to
an array of variable expressions.
"""
function gradient(O::Expression, vars::AbstractVector{<:Expression}; simplify = true)
    [expand_derivatives(Differential(v)(O),simplify) for v in vars]
end

"""
```julia
jacobian(ops::AbstractVector{<:Expression}, vars::AbstractVector{<:Expression}; simplify = true)
```

A helper function for computing the Jacobian of an array of expressions with respect to
an array of variable expressions.
"""
function jacobian(ops::AbstractVector{<:Expression}, vars::AbstractVector{<:Expression}; simplify = true)
    [expand_derivatives(Differential(v)(O),simplify) for O in ops, v in vars]
end

"""
```julia
hessian(O::Expression, vars::AbstractVector{<:Expression}; simplify = true)
```

A helper function for computing the Hessian of an expression with respect to
an array of variable expressions.
"""
function hessian(O::Expression, vars::AbstractVector{<:Expression}; simplify = true)
    [expand_derivatives(Differential(v2)(Differential(v1)(O)),simplify) for v1 in vars, v2 in vars]
end

function simplified_expr(O::Operation)
  if O isa Constant
    return O.value
  elseif isa(O.op, Differential)
    return :(derivative($(simplified_expr(O.args[1])),$(simplified_expr(O.op.x))))
  elseif isa(O.op, Variable)
    isempty(O.args) && return O.op.name
    return Expr(:call, Symbol(O.op), simplified_expr.(O.args)...)
  end
  if O.op === (^)
      if length(O.args) > 1  && O.args[2] isa Constant && O.args[2].value < 0
          return Expr(:call, :^, Expr(:call, :inv, simplified_expr(O.args[1])), -(O.args[2].value))
      end
  end
  return Expr(:call, Symbol(O.op), simplified_expr.(O.args)...)
end

simplified_expr(c::Constant) = c.value

function simplified_expr(eq::Equation)
    Expr(:(=), simplified_expr(eq.lhs), simplified_expr(eq.rhs))
end

simplified_expr(eq::AbstractArray) = simplified_expr.(eq)
