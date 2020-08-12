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

"""
```julia
jacobian_sparsity(ops::AbstractVector{<:Expression}, vars::AbstractVector{<:Expression})
```

Return the sparsity pattern of the Jacobian of an array of expressions with respect to
an array of variable expressions.
"""
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
    first_derivs = vec(jacobian([O], vars, simplify=simplify))
    n = length(vars)
    H = Array{Expression, 2}(undef,(n, n))
    fill!(H, 0)
    for i=1:n
        for j=1:i
            H[j, i] = H[i, j] = expand_derivatives(Differential(vars[i])(first_derivs[j]))
        end
    end
    H
end

isidx(x) = x isa TermCombination

"""
```julia
hessian_sparsity(ops::AbstractVector{<:Expression}, vars::AbstractVector{<:Expression})
```

Return the sparsity pattern of the Hessian of an array of expressions with respect to
an array of variable expressions.
"""
function hessian_sparsity end

let
    _scalar = one(TermCombination)

    linearity_propagator = [
          @rule +(~~xs) => reduce(+, filter(isidx, ~~xs), init=_scalar)
          @rule *(~~xs) => reduce(*, filter(isidx, ~~xs), init=_scalar)
          @rule (~f)(~x::(!isidx)) => _scalar
          @rule (~f)(~x::isidx) => if haslinearity_1(~f)
              combine_terms_1(linearity_1(~f), ~x)
          else
              error("Function of unknown linearity used: ", ~f)
          end
          @rule (^)(~x::isidx, ~y) => ~y isa Number && isone(~y) ? ~x : (~x) * (~x)
          @rule (~f)(~x, ~y) => begin
              if haslinearity_2(~f)
                  a = isidx(~x) ? ~x : _scalar
                  b = isidx(~y) ? ~y : _scalar
                  combine_terms_2(linearity(~f), a, b)
              else
                  error("Function of unknown linearity used: ", ~f)
              end
          end] |> Rewriters.Chain |> Rewriters.Postwalk |> Rewriters.Fixpoint

    global hessian_sparsity
    # we do this in a let block so that Revise works on the list of rules
    function hessian_sparsity(f, u)
        idx(i) = TermCombination(Set([Dict(i=>1)]))
        dict = Dict(SymbolicUtils.to_symbolic.(u) .=> idx.(1:length(u)))
        found = []
        f = Rewriters.Prewalk(x->haskey(dict, x) ? dict[x] : x)(to_symbolic(f))
        _sparse(linearity_propagator(to_symbolic(f)), length(u))
    end
end

"""
```julia
islinear(ex::Expression, u)
```
Check if an expression is linear with respect to a list of variable expressions.
"""
function islinear(ex::Expression, u)
    isempty(hessian_sparsity(ex, u).nzval)
end

"""
```julia
sparsehessian(O::Expression, vars::AbstractVector{<:Expression}; simplify = true)
```

A helper function for computing the sparse Hessian of an expression with respect to
an array of variable expressions.
"""
function sparsehessian(O::Expression, vars::AbstractVector{<:Expression}; simplify = true)
    S = hessian_sparsity(O, vars)
    I, J, _ = findnz(S)
    exprs = Array{Expression}(undef, length(I))
    fill!(exprs, 0)
    prev_j = 0
    d = nothing
    for (k, (i, j)) in enumerate(zip(I, J))
        j > i && continue
        if j != prev_j
            d = expand_derivatives(Differential(vars[j])(O), false)
        end
        expr = expand_derivatives(Differential(vars[i])(d), simplify)
        exprs[k] = expr
        prev_j = j
    end
    H = sparse(I, J, exprs, length(vars), length(vars))
    for (i, j) in zip(I, J)
        j > i && (H[i, j] = H[j, i])
    end
    return H
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
simplified_expr(x::Integer) = x
simplified_expr(x::AbstractFloat) = x
