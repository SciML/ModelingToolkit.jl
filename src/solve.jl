using SymbolicUtils: istree

function nterms(t)
    if istree(t)
        return reduce(+, map(nterms, arguments(t)), init=0)
    else
        return 1
    end
end
# Soft pivoted
# Note: we call this function with a matrix of Union{SymbolicUtils.Symbolic, Any}
function sym_lu(A)
    m, n = size(A)
    F = map(x->x isa Num ? x : Num(x), A)
    minmn = min(m, n)
    p = Vector{BlasInt}(undef, minmn)
    info = zero(BlasInt)
    for k = 1:minmn
        val, i = findmin(map(ii->_iszero(F[ii, k]) ? Inf : nterms(F[ii,k]), k:n))
        if !(val isa Symbolic) && (val isa Number) && val == Inf && iszero(info)
            info = k
        end
        i += k - 1
        # swap
        for j in 1:n
            F[k, j], F[i, j] = F[i, j], F[k, j]
        end
        p[k] = i
        for i in k+1:m
            F[i, k] = simplify(F[i, k] / F[k, k], polynorm=true)
        end
        for j = k+1:n
            for i in k+1:m
                F[i, j] = simplify(F[i, j] - F[i, k] * F[k, j], polynorm=true)
            end
        end
    end
    LU(F, p, info)
end

# Given a vector of equations and a
# list of independent variables,
# return the coefficient matrix `A` and a
# vector of constants (possibly symbolic) `b` such that
# A \ b will solve the equations for the vars
function A_b(eqs, vars)
    exprs = rhss(eqs) .- lhss(eqs)
    for ex in exprs
        @assert islinear(ex, vars)
    end
    A = jacobian(exprs, vars)
    b = A * vars - exprs
    A, b
end

"""
    solve_for(eqs::Vector, vars::Vector)

Solve the vector of equations `eqs` for a set of variables `vars`.

Assumes `length(eqs) == length(vars)`

Currently only works if all equations are linear.
"""
function solve_for(eqs, vars)
    A, b = A_b(eqs, vars)
    _solve(A, b)
end

function _solve(A, b)
    A = SymbolicUtils.simplify.(Num.(A), polynorm=true)
    b = SymbolicUtils.simplify.(Num.(b), polynorm=true)
    value.(SymbolicUtils.simplify.(sym_lu(A) \ b, polynorm=true))
end

# ldiv below

_iszero(x::Number) = iszero(x)
_isone(x::Number) = isone(x)
_iszero(::Term) = false
_isone(::Term) = false

function simplifying_dot(x,y)
    isempty(x) && return 0
    muls = map(x,y) do xi,yi
        _isone(xi) ? yi : _isone(yi) ? xi : _iszero(xi) ? 0 : _iszero(yi) ? 0 : xi * yi
    end

    reduce(muls) do acc, x
        _iszero(acc) ? x : _iszero(x) ? acc : acc + x
    end
end

function ldiv(F, b, x=b)
    L, U, p = F.L, F.U, F.p
    m, n = size(L)
    b = b[p]

    for i=n:-1:1
        sub = simplifying_dot(x[i+1:end], U[i,i+1:end])
        den = U[i,i]
        x[i] = _iszero(sub) ? b[i] : b[i] - sub
        x[i] = _isone(den) ? x[i] : _isone(-den) ? -x[i] : x[i] / den
    end

    # unit lower triangular solve first:
    for i=1:n
        sub = simplifying_dot(b[1:i-1], L[i, 1:i-1]) # this should be `b` not x
        x[i] = _iszero(sub) ? x[i] : x[i] - sub
    end
    x
end
