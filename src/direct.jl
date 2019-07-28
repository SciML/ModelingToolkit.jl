function gradient(O::Operation, vars::AbstractVector{Operation}; simplify = true)
    out = [expand_derivatives(Differential(v)(O)) for v in vars]
    simplify ? simplify_constants.(out) : out
end

function jacobian(ops::AbstractVector{Operation}, vars::AbstractVector{Operation}; simplify = true)
    out = [expand_derivatives(Differential(v)(O)) for O in ops, v in vars]
    simplify ? simplify_constants.(out) : out
end

function hessian(O::Operation, vars::AbstractVector{Operation}; simplify = true)
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

function simplified_expr(c::Constant)
    c.value
end

function simplified_expr(eq::Equation)
    Expr(:(=), simplified_expr(eq.lhs), simplified_expr(eq.rhs))
end

macro I(ex)
    name = :ICompile
    ret = return quote
        macro $(esc(name))()
            esc($ex)
        end
    end
end
