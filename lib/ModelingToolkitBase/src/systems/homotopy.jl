@register_symbolic homotopy(actual, simplified)

"""
    homotopy(actual, simplified)

Modelica homotopy operator (spec 3.7.4). At runtime and at L0 trivial rewrite
time, equivalent to `actual`. Future PRs may use `simplified` as a starting
point for parameter-sweep continuation during initialization.
"""
homotopy(actual::Real, simplified::Real) = actual

"""
    rewrite_trivial(expr)

Recursively replace every `homotopy(a, s)` node in `expr` with `a`. Preserves
`metadata` on every rebuilt node; `symtype` is inferred by `maketerm`.
Implemented as a hand-written walk over the TermInterface protocol — does NOT
use `@rule` (per AayushSabharwal feedback in #1385: per-node matcher overhead
is unacceptable at scale).
"""
function rewrite_trivial(expr)
    x = unwrap(expr)
    return _rewrite_trivial(x)
end

function _rewrite_trivial(x)
    if !iscall(x)
        return x
    end
    op = operation(x)
    args = arguments(x)
    new_args = map(_rewrite_trivial, args)
    if op === homotopy
        return new_args[1]
    end
    return maketerm(typeof(x), op, new_args, metadata(x))
end

"""
    has_homotopy(expr)

Return `true` iff `expr` contains at least one `homotopy(...)` node.
"""
function has_homotopy(expr)
    x = unwrap(expr)
    return _has_homotopy(x)
end

function _has_homotopy(x)
    if !iscall(x)
        return false
    end
    if operation(x) === homotopy
        return true
    end
    return any(_has_homotopy, arguments(x))
end

"""
    has_homotopy_in_equations(eqs)

Return `true` iff any equation in `eqs` (lhs or rhs) contains a `homotopy(...)`
node. `eqs` is an iterable of `Equation`.
"""
function has_homotopy_in_equations(eqs)
    for eq in eqs
        if has_homotopy(eq.lhs) || has_homotopy(eq.rhs)
            return true
        end
    end
    return false
end

"""
    rewrite_trivial_in_equations(eqs)

Return a new vector of `Equation`s with every `homotopy(a, s)` replaced by `a`
on both lhs and rhs. Original `eqs` not mutated; the system caller is
responsible for swapping the equation vector into the system.
"""
function rewrite_trivial_in_equations(eqs)
    return [Equation(rewrite_trivial(eq.lhs), rewrite_trivial(eq.rhs))
            for eq in eqs]
end
