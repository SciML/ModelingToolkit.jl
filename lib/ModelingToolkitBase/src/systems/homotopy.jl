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

"""
    rewrite_with_lambda(expr, λ = nothing)

Recursively replace every `homotopy(a, s)` node in `expr` with `(1 - λ)*s + λ*a`,
where `λ` is a single shared parameter for the whole expression (allocated lazily
if not supplied, default name `__homotopy_λ`, default value `1.0`).
Returns `(new_expr, λ)`.

At λ=1 the lowered expression reduces numerically to `actual` (trivial form);
at λ=0 it reduces to `simplified`. `HomotopySweep` walks λ from 0 → 1.

Hand-written TermInterface walk; no `@rule` (per AayushSabharwal #1385 feedback).
"""
function rewrite_with_lambda(expr, λ = nothing)
    if λ === nothing
        λ = only(@parameters __homotopy_λ = 1.0)
    end
    x = unwrap(expr)
    return _rewrite_with_lambda(x, λ), λ
end

function _rewrite_with_lambda(x, λ)
    if !iscall(x)
        return x
    end
    op = operation(x)
    args = arguments(x)
    new_args = map(arg -> _rewrite_with_lambda(arg, λ), args)
    if op === homotopy
        a, s = new_args[1], new_args[2]
        return (1 - λ) * s + λ * a
    end
    return maketerm(typeof(x), op, new_args, metadata(x))
end

"""
    rewrite_with_lambda_in_equations(eqs, λ = nothing)

Return `(new_eqs, λ)`. All `homotopy(a, s)` nodes across all equations are
replaced using a single shared `λ` parameter. Original `eqs` not mutated.
"""
function rewrite_with_lambda_in_equations(eqs, λ = nothing)
    if λ === nothing
        λ = only(@parameters __homotopy_λ = 1.0)
    end
    new_eqs = [Equation(_rewrite_with_lambda(unwrap(eq.lhs), λ),
                        _rewrite_with_lambda(unwrap(eq.rhs), λ)) for eq in eqs]
    return new_eqs, λ
end

"""
    has_any_homotopy(sys)

Return `true` iff `sys` contains a `homotopy(...)` node anywhere in its
equations OR its observed equations. Used by `add_homotopy_parameter` to
decide whether to inject `__homotopy_λ` into a parent system at `complete`
time.
"""
function has_any_homotopy(sys)
    has_homotopy_in_equations(equations(sys)) && return true
    obs = observed(sys)
    obs === nothing && return false
    for eq in obs
        (has_homotopy(eq.lhs) || has_homotopy(eq.rhs)) && return true
    end
    return false
end

"""
    add_homotopy_parameter(sys)

If `sys` contains `homotopy(a, s)` nodes (in equations or observed), lower
every node to `(1 - λ)*s + λ*a` and inject a shared `__homotopy_λ` parameter
(default value `1.0`) into the system's parameter list. Mirrors the
`add_initialization_parameters` pattern (`abstractsystem.jl:559`) — same
lifecycle point inside `complete()`, so the parent and the init system end
up sharing one identity for `__homotopy_λ`. This lets downstream codepaths
like `MTKParametersReconstructor` resolve λ when it appears in init-system
observed expressions.

At λ=1 (the default) the lowered expression reduces numerically to `actual`;
`HomotopySweep` walks λ from 0 → 1 to obtain the actual root from a
`simplified` starting point.

No-op if `sys` doesn't support initialization (Modelica spec 3.7.4 only
applies during init), is itself an initializesystem, or contains no
`homotopy(...)` nodes.
"""
function add_homotopy_parameter(sys::AbstractSystem)
    supports_initialization(sys) || return sys
    is_initializesystem(sys) && return sys
    has_any_homotopy(sys) || return sys

    λ = only(@parameters __homotopy_λ = 1.0)
    λ_sym = unwrap(λ)

    # Lower homotopy nodes in equations
    new_eqs, _ = rewrite_with_lambda_in_equations(equations(sys), λ_sym)
    @set! sys.eqs = new_eqs

    # Lower homotopy nodes in observed equations (PressureDrop-style: m_flow
    # is eliminated by `mtkcompile` and its definition lives in observed —
    # downstream observed codegen must see lowered form, not opaque homotopy).
    obs = observed(sys)
    if obs !== nothing && !isempty(obs)
        new_obs = [Equation(_rewrite_with_lambda(unwrap(eq.lhs), λ_sym),
                            _rewrite_with_lambda(unwrap(eq.rhs), λ_sym)) for eq in obs]
        @set! sys.observed = new_obs
    end

    # Inject λ as a parameter on the parent — idempotent
    current_ps = get_ps(sys)
    if !any(p -> isequal(p, λ_sym), current_ps)
        @set! sys.ps = vcat(current_ps, λ_sym)
    end

    return sys
end
