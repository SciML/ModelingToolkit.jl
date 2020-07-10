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
sparsejacobian(ops::AbstractVector{<:Expression}, vars::AbstractVector{<:Expression}; simplify = true)
```

A helper function for computing the sparse Jacobian of an array of expressions with respect to
an array of variable expressions.
"""
function sparsejacobian(ops::AbstractVector{<:Expression}, vars::AbstractVector{<:Expression}; simplify = true)
    I = Int[]
    J = Int[]
    du = Expression[]

    sp = jacobian_sparsity(ops, vars)
    I,J,_ = findnz(sp)

    exprs = Expression[]

    for (i,j) in zip(I, J)
        push!(exprs, expand_derivatives(Differential(vars[j])(ops[i]), simplify))
    end
    sparse(I,J, exprs, length(ops), length(vars))
end

function jacobian_sparsity(du, u)
    dict = Dict(zip(to_symbolic.(u), 1:length(u)))

    i = Ref(1)
    I = Int[]
    J = Int[]

    # This rewriter notes down which u's appear in a
    # given du (whose index is stored in the `i` Ref)
    r = [@rule ~x::(x->haskey(dict, x)) => begin
        push!(I, i[])
        push!(J, dict[~x])
        nothing
    end] |> Rewriters.Chain |> Rewriters.Postwalk

    for ii = 1:length(du)
        i[] = ii
        r(to_symbolic(du[ii]))
    end

    sparse(I, J, true, length(du), length(u))
end

"""
```julia
hessian(O::Expression, vars::AbstractVector{<:Expression}; simplify = true)
```

A helper function for computing the Hessian of an expression with respect to
an array of variable expressions.
"""
function hessian(O::Expression, vars::AbstractVector{<:Expression}; simplify = true)
    jacobian(vec(jacobian([O], vars, simplify=simplify)), vars, simplify=simplify)
end

isidx(x) = x isa TermCombination

"""
```julia
sparsehessian(O::Expression, vars::AbstractVector{<:Expression}; simplify = true)
```

A helper function for computing the sparse Hessian of an expression with respect to
an array of variable expressions.
"""
function hessian_sparsity(f, u)
    idx(i) = TermCombination(Set([Dict(i=>1)]))
    dict = Dict(SymbolicUtils.to_symbolic.(u) .=> idx.(1:length(u)))
    found = []
    f = Rewriters.Prewalk(Rewriters.Chain([@rule ~x::(x->haskey(dict, x)) => dict[~x]]))(to_symbolic(f))

    # condense
    z = one(TermCombination)
    rr = [@rule +(~~xs) => reduce(+, filter(isidx, ~~xs), init=z)
          @rule *(~~xs) => reduce(*, filter(isidx, ~~xs), init=z)
          @rule (~f)(~x::(!isidx)) => z
          @rule (~f)(~x::isidx) => if haslinearity(~f, Val{1}())
              combine_terms(linearity(~f, Val{1}()), ~x)
          else
              error("Function of unknown linearity used ", ~f)
          end
          @rule (~f)(~x, ~y) => begin
              if haslinearity(~f, Val{2}())
                  a = isidx(~x) ? ~x : z
                  b = isidx(~y) ? ~y : z
                  combine_terms(linearity(~f, Val{2}()), a, b)
              else
                  error("Function of unknown linearity used ", ~f)
              end
          end]

    _sparse(Rewriters.Fixpoint(Rewriters.Postwalk(Rewriters.Chain(rr)))(to_symbolic(f)), length(u))
end

function sparsehessian(O::Expression, vars::AbstractVector{<:Expression}; simplify = true)
    S = hessian_sparsity(O, vars)
    I, J, _ = findnz(S)
    exprs = Expression[]
    prev_j = 0
    d = nothing
    for (i, j) in zip(I, J)
        if j != prev_j
            d = expand_derivatives(Differential(vars[j])(O), false)
        end
        expr = expand_derivatives(Differential(vars[i])(d), simplify)
        push!(exprs, expr)
        prev_j = j
    end
    sparse(I, J, exprs, length(vars), length(vars))
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
