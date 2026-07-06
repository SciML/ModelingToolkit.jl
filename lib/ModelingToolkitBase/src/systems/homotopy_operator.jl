# Modelica's `homotopy(actual, simplified)` operator (spec 3.7.4.2) and its
# lowering to a continuation residual. NOT to be confused with
# `systems/nonlinear/homotopy_continuation.jl`, which is the unrelated polynomial
# HomotopyContinuation.jl rooting feature.
@register_symbolic homotopy(actual, simplified)

"""
    homotopy(actual, simplified)

The Modelica homotopy operator (Modelica Specification 3.7.4.2). Annotating an
expression as `homotopy(actual, simplified)` declares that `simplified` is an
easy-to-solve approximation of `actual`: a continuation (homotopy) solver can
start from the `simplified` equations and continuously deform them into the
`actual` ones. This stabilises the solution of nonlinear systems that are hard
to solve from a cold start — pressure-driven flow networks, power-flow
equations, chemical equilibria, equations with multiple roots — by walking from
the easy solution to the true one. Wherever the operator is evaluated
numerically (outside a continuation solve) it is simply `actual`, as the
Modelica spec prescribes; initialization is one common consumer, not the only
one.

# Arguments

  - `actual`: the real expression. Used wherever the operator is evaluated
    numerically, and reached at the end of the continuation (`λ = 1`).
  - `simplified`: an approximation of `actual` whose solution is easy to reach
    from the available guess. Only influences the continuation path.

Both arguments are scalar `Real` expressions, matching Modelica's restriction of
`homotopy` to scalar `Real`s; the runtime method is
`homotopy(actual::Real, simplified::Real) = actual`. For array equations,
broadcast the operator elementwise — `homotopy.(actual, simplified)` — which
creates one `homotopy` node per element; the continuation lowering rewrites each
node independently, and all of them share the single continuation parameter `λ`.

# Behavior

`homotopy(actual, simplified)` stays an opaque symbolic operator through `System`
construction, `mtkcompile`, and runtime code generation — no continuation
parameter is injected into the system, and systems that do not use the operator
go through a byte-identical pipeline. Symbolic differentiation works through the
operator: nodewise derivative rules keep symbolic jacobians, `tgrad`, and index
reduction consistent — at runtime differentiated equations reproduce `actual`'s
derivative, and along the continuation they follow the blended expression below.

Building a [`SciMLBase.HomotopyProblem`](@ref) from a system that contains
`homotopy` nodes — directly with `HomotopyProblem(sys, op)`, or automatically
through `AbstractNonlinearProblem(sys, op)` — regenerates the residual with every
`homotopy(actual, simplified)` replaced by the convex blend

```julia
(1 - λ) * simplified + λ * actual
```

compiled as `f(u, p, λ)` with `λ` an explicit trailing argument — `λ` is never
added to the system's parameters, and the user's parameter object `p` passes
through untouched. All `homotopy` calls in a system share the single `λ`, per the
Modelica spec's recommendation of (conceptually) one homotopy iteration over the
whole model. The resulting `HomotopyProblem` (with `λspan` defaulting to
`(0.0, 1.0)`) can be solved with any algorithm that supports it; `solve(prob)`
with no algorithm picks a default that sweeps `λ` from `0` (`simplified`) to
`1` (`actual`).

# Example

The equation `0 = atan(y - 3)` has its root at `y = 3`, but Newton from `y = 12`
diverges because `atan` saturates. Annotating with `simplified = y` (root at
`y = 0`) lets the continuation walk to the true root:

```julia
using ModelingToolkit, NonlinearSolve

@variables y
@mtkcompile sys = System([0 ~ homotopy(atan(y - 3), y)])
prob = HomotopyProblem(sys, [y => 12.0])
sol = solve(prob)
sol[y] # ≈ 3.0 — the continuation rescued the out-of-basin guess
```

# Notes

  - Runtime cost: outside a continuation solve the generated code calls the
    numeric fallback, so the `simplified` argument expression is evaluated and
    its value discarded (arguments are evaluated before the call). This small
    overhead is borne only by systems that use `homotopy`; all other systems are
    unaffected.
  - See the [Homotopy](@ref homotopy) documentation page for construction details
    and the continuation solver's options.

**Reference:** Modelica Specification 3.7.4.2
<https://specification.modelica.org/master/operators-and-expressions.html#homotopy>
"""
homotopy(actual::Real, simplified::Real) = actual

# Nodewise symbolic derivatives (Modelica-faithful).
#
# ∂homotopy(a, s)/∂a = homotopy(1, 0) and ∂homotopy(a, s)/∂s = homotopy(0, 1),
# as OPAQUE nodes — deliberately NOT folded to the constants 1/0. The chain
# rule then yields `homotopy(1, 0)*a′ + homotopy(0, 1)*s′`, so that:
#   * at runtime the numeric fallback folds the factors to 1/0, reproducing
#     exactly the derivative of `actual` (Modelica spec: homotopy ≡ actual
#     outside a continuation solve);
#   * the continuation rewrite `homotopy(a, s) → (1-λ)s + λa` blends the same
#     factors to λ/(1-λ), i.e. the true derivative of the blended expression —
#     this keeps symbolic jacobians and index reduction (dummy derivatives)
#     consistent with OMC's homotopy lowering.
# `@register_derivative` requires the rule to return a symbolic result;
# `SymbolicUtils.term` builds the opaque call node without invoking the numeric
# fallback on the constant arguments, and the constants are wrapped in
# `Symbolics.SConst` to construct the node explicitly from symbolic constants.
Symbolics.@register_derivative homotopy(actual, simplified) 1 SymbolicUtils.term(
    homotopy, Symbolics.SConst(1), Symbolics.SConst(0); type = Real
)
Symbolics.@register_derivative homotopy(actual, simplified) 2 SymbolicUtils.term(
    homotopy, Symbolics.SConst(0), Symbolics.SConst(1); type = Real
)

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
    has_any_homotopy(sys)

Return `true` iff `sys` contains a `homotopy(...)` node anywhere in its
equations OR its observed equations.
"""
function has_any_homotopy(sys)
    has_homotopy_in_equations(equations(sys)) && return true
    # `observed(sys)` always returns a `Vector{Equation}`, never `nothing`.
    return any(eq -> has_homotopy(eq.lhs) || has_homotopy(eq.rhs), observed(sys))
end

"""
    _rewrite_with_lambda(x, λ)

Recursively replace every `homotopy(a, s)` node in the unwrapped expression `x`
with `(1 - λ)*s + λ*a`, where `λ` is supplied by the caller. At λ=1 the lowered
expression reduces numerically to `actual` (trivial form); at λ=0 it reduces to
`simplified`.

Hand-written TermInterface walk rather than a SymbolicUtils `@rule` — keeps the
recursion explicit and avoids the rule-rewriter overhead for a single-node rewrite.
"""
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

# The continuation parameter is a fixed sentinel symbol, NOT `gensym`. A
# `gensym` name embeds a process-global counter, so the same system would lower
# to a different name (and thus a different generated `Expr`) in the precompile
# process than in the user session — defeating the `RuntimeGeneratedFunctions`
# Expr-hash cache and breaking precompilation. A fixed name makes codegen a pure
# function of the system. The `ₘₜₖ` unicode suffix (matching MTK's other
# internal sentinels, e.g. `__log_assertions_ₘₜₖ`) makes a collision with a
# user-chosen symbol practically impossible; `lower_homotopy` additionally
# guards against it, since `Sym` equality is name-keyed and a collision would
# otherwise silently rewrite `λ` into a state/parameter access and drop the
# `λ`-dependence of the residual.
const HOMOTOPY_LAMBDA = unwrap(only(@variables __homotopy_λₘₜₖ))

"""
    lower_homotopy(sys::AbstractSystem) -> (shadow, λ)

Independent lowering pass that feeds homotopy-swept residual code generation.
Returns a shadow copy of `sys` in which every `homotopy(actual, simplified)`
node in the equations and observed equations is replaced by the convex blend
`(1 - λ)*simplified + λ*actual`, plus the single shared bare symbolic `λ`
([`HOMOTOPY_LAMBDA`](@ref)).

`sys` must be a `complete`d (flattened) system: the pass reads equations and
observed through accessors but reconstructs the toplevel fields, which is only
coherent once subsystems have been flattened away.

`λ` is a fixed reserved sentinel and is **not** added to the system's parameters
— it is intended to become an explicit trailing argument of the swept residual
rather than an entry of `p`, per `SciMLBase.HomotopyProblem`'s `f(u, p, λ)`
contract. The shadow never passes through `mtkcompile` (`TearingState` would
reject the free λ); it exists only to feed the swept-residual code generation.
"""
function lower_homotopy(sys::AbstractSystem)
    if !iscomplete(sys)
        throw(
            ArgumentError(
                "`lower_homotopy` expects a `complete`d (flattened) system."
            )
        )
    end
    λ = HOMOTOPY_LAMBDA
    λname = getname(λ)
    # Guard against a user symbol named exactly like the sentinel: a collision in
    # the unknowns/parameters would silently rewrite λ into a state/parameter
    # access and drop the residual's λ-dependence. Observed lhs are included for
    # completeness (a variable eliminated into observed by `mtkcompile`).
    for v in Iterators.flatten(
            (
                unknowns(sys), parameters(sys; initial_parameters = true),
                (eq.lhs for eq in observed(sys)),
            )
        )
        if hasname(v) && getname(v) === λname
            throw(
                ArgumentError(
                    "the homotopy continuation sentinel `$(λname)` collides with a " *
                        "symbol of the system; this name is reserved for internal use."
                )
            )
        end
    end
    rew(x) = _rewrite_with_lambda(unwrap(x), λ)
    new_eqs = [Equation(rew(eq.lhs), rew(eq.rhs)) for eq in equations(sys)]
    new_obs = [Equation(rew(eq.lhs), rew(eq.rhs)) for eq in observed(sys)]
    # All other fields (incl. `tearing_state`/`schedule`) are carried over by
    # reference and remain valid: the rewrite is structure-preserving — same
    # equation count/order/incidence, modulo the free λ.
    return ConstructionBase.setproperties(sys; eqs = new_eqs, observed = new_obs), λ
end

"""
    generate_homotopy_residual(shadow, λ; kwargs...)

Compile the swept residual `f(u, p, λ)` (out-of-place) / `f(du, u, p, λ)`
(in-place) from the homotopy-lowered shadow system. `λ` is the bare symbol
returned by [`lower_homotopy`](@ref); it is appended as an explicit trailing
scalar argument (the "t slot" of the 4-argument convention used by
`SciMLBase.HomotopyProblem`), and the user's parameter object is passed through
untouched (λ is excluded from the `MTKParameters` destructuring).

This is exactly [`generate_rhs`](@ref) for the time-independent case with `λ`
threaded through its `extra_args`, behind two loud guards: both misuse modes
(a non-lowered system with raw `homotopy(...)` nodes, or a time-dependent
system) were probed to fail SILENTLY with wrong numerics rather than erroring.

# Keyword Arguments

$GENERATE_X_KWARGS

All other keyword arguments are forwarded to [`generate_rhs`](@ref).
"""
function generate_homotopy_residual(
        shadow::AbstractSystem, λ; expression = Val{false}, wrap_gfw = Val{true},
        kwargs...
    )
    if has_any_homotopy(shadow)
        throw(
            ArgumentError(
                "`generate_homotopy_residual` expects a homotopy-lowered shadow " *
                    "system as returned by `lower_homotopy`; this system still " *
                    "contains raw `homotopy(...)` nodes, which codegen would fold " *
                    "to `actual` via the numeric fallback, producing a " *
                    "λ-independent residual."
            )
        )
    end
    if is_time_dependent(shadow)
        throw(
            ArgumentError(
                "`generate_homotopy_residual` expects a time-independent system; " *
                    "got a time-dependent one."
            )
        )
    end
    return generate_rhs(
        shadow; extra_args = (λ,), expression, wrap_gfw, kwargs...
    )
end
