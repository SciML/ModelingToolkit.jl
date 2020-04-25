"""
```julia
gradient(O::Expression, vars::AbstractVector{<:Expression}; simplify = true)
```

A helper function for computing the gradient of an expression with respect to
an array of variable expressions.
"""
function gradient(O::Expression, vars::AbstractVector{<:Expression}; simplify = true)
    out = [expand_derivatives(Differential(v)(O)) for v in vars]
    simplify ? simplify_constants.(out) : out
end

"""
```julia
jacobian(O::Expression, vars::AbstractVector{<:Expression}; simplify = true)
```

A helper function for computing the Jacobian of an array of expressions with respect to
an array of variable expressions.
"""
function jacobian(ops::AbstractVector{<:Expression}, vars::AbstractVector{<:Expression}; simplify = true)
    out = [expand_derivatives(Differential(v)(O)) for O in ops, v in vars]
    simplify ? simplify_constants.(out) : out
end

"""
```julia
hessian(O::Expression, vars::AbstractVector{<:Expression}; simplify = true)
```

A helper function for computing the Hessian of an expression with respect to
an array of variable expressions.
"""
function hessian(O::Expression, vars::AbstractVector{<:Expression}; simplify = true)
    out = [expand_derivatives(Differential(v2)(Differential(v1)(O))) for v1 in vars, v2 in vars]
    simplify ? simplify_constants.(out) : out
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
  return Expr(:call, Symbol(O.op), simplified_expr.(O.args)...)
end

simplified_expr(c::Constant) = c.value
simplified_expr(c) = c

function simplified_expr(eq::Equation)
    Expr(:(=), simplified_expr(eq.lhs), simplified_expr(eq.rhs))
end

simplified_expr(eq::AbstractArray) = simplified_expr.(eq)
