"""
```julia
derivative(O, v; simplify = true)
```

A helper function for computing the derivative of an expression with respect to
`var`.
"""
function derivative(O, v; simplify = true)
    if O isa AbstractArray
        Num[Num(expand_derivatives(Differential(v)(value(o)), simplify)) for o in O]
    else
        Num(expand_derivatives(Differential(v)(value(O)), simplify))
    end
end

"""
```julia
gradient(O, vars::AbstractVector; simplify = true)
```

A helper function for computing the gradient of an expression with respect to
an array of variable expressions.
"""
function gradient(O, vars::AbstractVector; simplify = true)
    Num[Num(expand_derivatives(Differential(v)(value(O)),simplify)) for v in vars]
end

"""
```julia
jacobian(ops::AbstractVector, vars::AbstractVector; simplify = true)
```

A helper function for computing the Jacobian of an array of expressions with respect to
an array of variable expressions.
"""
function jacobian(ops::AbstractVector, vars::AbstractVector; simplify = true)
    Num[Num(expand_derivatives(Differential(value(v))(value(O)),simplify)) for O in ops, v in vars]
end

"""
```julia
sparsejacobian(ops::AbstractVector, vars::AbstractVector; simplify = true)
```

A helper function for computing the sparse Jacobian of an array of expressions with respect to
an array of variable expressions.
"""
function sparsejacobian(ops::AbstractVector, vars::AbstractVector; simplify = true)
    I = Int[]
    J = Int[]
    du = Num[]

    sp = jacobian_sparsity(ops, vars)
    I,J,_ = findnz(sp)

    exprs = Num[]

    for (i,j) in zip(I, J)
        push!(exprs, Num(expand_derivatives(Differential(vars[j])(ops[i]), simplify)))
    end
    sparse(I,J, exprs, length(ops), length(vars))
end

"""
```julia
jacobian_sparsity(ops::AbstractVector, vars::AbstractVector)
```

Return the sparsity pattern of the Jacobian of an array of expressions with respect to
an array of variable expressions.
"""
function jacobian_sparsity(du, u)
    du = map(value, du)
    u = map(value, u)
    dict = Dict(zip(u, 1:length(u)))

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
        r(du[ii])
    end

    sparse(I, J, true, length(du), length(u))
end

"""
    exprs_occur_in(exprs::Vector, expr)

Return an array of booleans `finds` where `finds[i]` is true if `exprs[i]` occurs in `expr`
false otherwise.
"""
function exprs_occur_in(exprs, expr)
    vec(jacobian_sparsity([expr], exprs))
end

"""
```julia
hessian(O, vars::AbstractVector; simplify = true)
```

A helper function for computing the Hessian of an expression with respect to
an array of variable expressions.
"""
function hessian(O, vars::AbstractVector; simplify = true)
    vars = map(value, vars)
    first_derivs = map(value, vec(jacobian([values(O)], vars, simplify=simplify)))
    n = length(vars)
    H = Array{Num, 2}(undef,(n, n))
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
hessian_sparsity(ops::AbstractVector, vars::AbstractVector)
```

Return the sparsity pattern of the Hessian of an array of expressions with respect to
an array of variable expressions.
"""
function hessian_sparsity end

let
    # we do this in a let block so that Revise works on the list of rules

    _scalar = one(TermCombination)

    simterm(t, f, args) = Term{Any}(f, args)
    linearity_rules = [
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
                  combine_terms_2(linearity_2(~f), a, b)
              else
                  error("Function of unknown linearity used: ", ~f)
              end
          end]
    linearity_propagator = Fixpoint(Postwalk(Chain(linearity_rules); similarterm=simterm))

    global hessian_sparsity

    """
    ```julia
    hessian_sparsity(ops::AbstractVector, vars::AbstractVector)
    ```

    Return the sparsity pattern of the Hessian of an array of expressions with respect to
    an array of variable expressions.
    """
    function hessian_sparsity(f, u)
        @assert !(f isa AbstractArray)
        f = value(f)
        u = map(value, u)
        idx(i) = TermCombination(Set([Dict(i=>1)]))
        dict = Dict(u .=> idx.(1:length(u)))
        f = Rewriters.Prewalk(x->haskey(dict, x) ? dict[x] : x; similarterm=simterm)(f)
        lp = linearity_propagator(f)
        _sparse(lp, length(u))
    end
end

"""
```julia
islinear(ex, u)
```
Check if an expression is linear with respect to a list of variable expressions.
"""
function islinear(ex, u)
    isempty(hessian_sparsity(ex, u).nzval)
end

"""
```julia
sparsehessian(O, vars::AbstractVector; simplify = true)
```

A helper function for computing the sparse Hessian of an expression with respect to
an array of variable expressions.
"""
function sparsehessian(O, vars::AbstractVector; simplify = true)
    O = value(O)
    vars = map(value, vars)
    S = hessian_sparsity(O, vars)
    I, J, _ = findnz(S)
    exprs = Array{Num}(undef, length(I))
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

"""
    toexpr(O::Union{Symbolics,Num,Equation,AbstractArray}; canonicalize=true) -> Expr

Convert `Symbolics` into `Expr`. If `canonicalize`, then we turn exprs like
`x^(-n)` into `inv(x)^n` to avoid type error when evaluating.
"""
function toexpr(O; canonicalize=true)
    if canonicalize
        canonical, O = canonicalexpr(O)
        canonical && return O
    else
        !istree(O) && return O
    end

    op = operation(O)
    args = arguments(O)
    if op isa Differential
        ex = toexpr(args[1]; canonicalize=canonicalize)
        wrt = toexpr(op.x; canonicalize=canonicalize)
        return :(Differential($wrt)($ex))
    elseif op isa Sym
        isempty(args) && return nameof(op)
        return Expr(:call, toexpr(op; canonicalize=canonicalize), toexpr(args; canonicalize=canonicalize)...)
    end
    return Expr(:call, op, toexpr(args; canonicalize=canonicalize)...)
end
toexpr(s::Sym; kw...) = nameof(s)

"""
    canonicalexpr(O) -> (canonical::Bool, expr)

Canonicalize `O`. Return `canonical` if `expr` is valid code to generate.
"""
function canonicalexpr(O)
    !istree(O) && return true, O
    op = operation(O)
    args = arguments(O)
    if op === (^)
        if length(args) == 2 && args[2] isa Number && args[2] < 0
            ex = toexpr(args[1])
            if args[2] == -1
                expr = Expr(:call, inv, ex)
            else
                expr = Expr(:call, ^, Expr(:call, inv, ex), -args[2])
            end
            return true, expr
        end
    end
    return false, O
end

function toexpr(eq::Equation; kw...)
    Expr(:(=), toexpr(eq.lhs; kw...), toexpr(eq.rhs; kw...))
end

toexpr(eqs::AbstractArray; kw...) = map(eq->toexpr(eq; kw...), eqs)
toexpr(x::Integer; kw...) = x
toexpr(x::AbstractFloat; kw...) = x
toexpr(x::Num; kw...) = toexpr(value(x); kw...)
