# Modelica homotopy operator — L0 trivial form (spec 3.7.4).
# `homotopy(actual, simplified)` is a symbolic primitive; at runtime and in
# the init pipeline it evaluates to `actual`. The `simplified` argument is
# carried symbolically so future PRs (L1 parameter sweep) can use it for
# continuation without changing call sites.

using Symbolics
using SymbolicUtils

@register_symbolic homotopy(actual, simplified)

"""
    homotopy(actual, simplified)

Modelica homotopy operator (spec 3.7.4). At runtime and at L0 trivial rewrite
time, equivalent to `actual`. Future PRs may use `simplified` as a starting
point for parameter-sweep continuation during initialization.
"""
homotopy(actual::Real, simplified::Real) = actual
homotopy(actual::AbstractArray, simplified::AbstractArray) = actual

"""
    rewrite_trivial(expr)

Recursively replace every `homotopy(a, s)` node in `expr` with `a`. Preserves
`metadata` and `symtype` on every rebuilt node. Implemented as a hand-written
walk over `TermInterface` — does NOT use `@rule` (per AayushSabharwal feedback
in #1385: per-node matcher overhead is unacceptable at scale).
"""
function rewrite_trivial(expr)
    x = Symbolics.unwrap(expr)
    return _rewrite_trivial(x)
end

function _rewrite_trivial(x)
    if !SymbolicUtils.iscall(x)
        return x
    end
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    new_args = map(_rewrite_trivial, args)
    if op === homotopy
        length(new_args) == 2 || throw(ArgumentError(
            "homotopy() expects exactly 2 arguments, got $(length(new_args))"))
        return new_args[1]
    end
    return SymbolicUtils.maketerm(
        typeof(x), op, new_args,
        SymbolicUtils.metadata(x),
    )
end

"""
    has_homotopy(expr)

Return `true` iff `expr` contains at least one `homotopy(...)` node.
"""
function has_homotopy(expr)
    x = Symbolics.unwrap(expr)
    return _has_homotopy(x)
end

function _has_homotopy(x)
    if !SymbolicUtils.iscall(x)
        return false
    end
    if SymbolicUtils.operation(x) === homotopy
        return true
    end
    return any(_has_homotopy, SymbolicUtils.arguments(x))
end
