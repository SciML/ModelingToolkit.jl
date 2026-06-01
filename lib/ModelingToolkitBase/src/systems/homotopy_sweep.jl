using CommonSolve: CommonSolve
using SciMLBase: SciMLBase, NonlinearProblem, remake,
                 successful_retcode, state_values, parameter_values
using Setfield: @set
using SimpleNonlinearSolve: SimpleNewtonRaphson

# An inner residual/continuation solve must never re-trigger DAE initialization.
# When the parent is built with an explicit `initializealg`, that `OverrideInit`
# is baked into the init problem's own `prob.kwargs` and `solve_up` merges it back
# even if we drop it from the forwarded call kwargs. NonlinearSolveBase then runs
# its `OverrideInit` init path on the inner problem — whose `f.initialization_data`
# is `nothing` — and throws a `FieldError` (that path is unguarded, unlike its
# default init path). Force `initializealg = NoInit()` on the inner solve so the
# call kwarg overrides any baked-in `OverrideInit` and no inner re-init runs.
_inner_solve_kwargs(kwargs) = (;
    Base.structdiff((; kwargs...), NamedTuple{(:initializealg,)})...,
    initializealg = SciMLBase.NoInit())

"""
    HomotopySweep(; inner = NewtonRaphson(), schedule = 0.0:0.1:1.0,
                  set_λ!, maxiters_per_step = nothing)

Parameter-sweep continuation algorithm for `NonlinearProblem` lowered by
`rewrite_with_lambda`. Each step `λ_k ∈ schedule` is solved by `inner` with
`u0` warm-started from the previous step's solution. λ is written into the
problem's parameter vector via `set_λ!`, normally a `SymbolicIndexingInterface.setp`
handle (or a `(p, v) -> new_p` closure in tests).

Failure policy: if any step's inner solve returns an unsuccessful retcode,
`HomotopySweep` returns that step's solution unchanged. No auto-fallback.
Auto-fallback to trivial is the job of `TrivialThenSweep`, not this struct.
"""
struct HomotopySweep{Inner, Sched, Setter} <: SciMLBase.AbstractNonlinearAlgorithm
    inner::Inner
    schedule::Sched
    set_λ!::Setter
    maxiters_per_step::Union{Int, Nothing}
end

function HomotopySweep(; inner = nothing, schedule = 0.0:0.1:1.0,
                        set_λ!, maxiters_per_step = nothing)
    return HomotopySweep(inner === nothing ? _default_inner() : inner,
                         schedule, set_λ!, maxiters_per_step)
end

# Default inner solver for the homotopy continuation algorithms. NonlinearSolve
# is NOT a runtime `[deps]` of ModelingToolkitBase, so its richer `NewtonRaphson`
# cannot be referenced directly here. The `MTKNonlinearSolveExt` package
# extension (a `[weakdeps]` extension) populates this Ref with a thunk returning
# `NonlinearSolve.NewtonRaphson()` whenever NonlinearSolve is loaded. Until then
# `_default_inner` falls back to `SimpleNonlinearSolve.SimpleNewtonRaphson` (a
# runtime `[deps]`) so a system carrying `homotopy(...)` can be built and solved
# out of the box. Load NonlinearSolve (`using NonlinearSolve`) to restore its
# `NewtonRaphson` as the default, or pass `inner = ...` explicitly.
const _DEFAULT_INNER_FACTORY = Ref{Union{Nothing, Function}}(nothing)

function _default_inner()
    factory = _DEFAULT_INNER_FACTORY[]
    return factory === nothing ? SimpleNewtonRaphson() : factory()
end

function CommonSolve.solve(prob::NonlinearProblem, alg::HomotopySweep; kwargs...)
    u_curr = copy(state_values(prob))
    last_sol = nothing
    # `set_λ!` is expected to be a `SymbolicIndexingInterface.setp` handle in
    # production (mutates `prob.p` in place via `setindex!`, returns nothing).
    # Tests pass an OOP `(p, v) -> new_p` closure; both paths are supported
    # by the `ret === nothing ? p_in : ret` normalization below.
    for λ in alg.schedule
        p_in = copy(parameter_values(prob))
        ret = alg.set_λ!(p_in, λ)
        new_p = ret === nothing ? p_in : ret
        step_prob = remake(prob; u0 = u_curr, p = new_p)
        base_kwargs = _inner_solve_kwargs(kwargs)
        inner_kwargs = alg.maxiters_per_step === nothing ?
                       base_kwargs : (; base_kwargs..., maxiters = alg.maxiters_per_step)
        step_sol = CommonSolve.solve(step_prob, alg.inner; inner_kwargs...)
        last_sol = step_sol
        if !successful_retcode(step_sol)
            return step_sol
        end
        u_curr = copy(step_sol.u)
    end
    return last_sol
end

"""
    TrivialHomotopy(; inner = NewtonRaphson())

Trivial-form init: leaves `__homotopy_λ` at its default 1.0 (so the lowered
system is `actual`) and solves once with `inner`. Equivalent to OMC's
`-noHomotopyOnFirstTry`. Use when the guess is reliably close to the actual
root and sweep cost is unwarranted.
"""
struct TrivialHomotopy{Inner} <: SciMLBase.AbstractNonlinearAlgorithm
    inner::Inner
end

TrivialHomotopy(; inner = nothing) =
    TrivialHomotopy(inner === nothing ? _default_inner() : inner)

function CommonSolve.solve(prob::NonlinearProblem, alg::TrivialHomotopy; kwargs...)
    return CommonSolve.solve(prob, alg.inner; _inner_solve_kwargs(kwargs)...)
end

"""
    TrivialThenSweep(; trivial = TrivialHomotopy(), sweep = HomotopySweep(...))

Composite algorithm matching OpenModelica's default user experience: attempt
the trivial (single-Newton) solve first, and on unsuccessful retcode fall
back to a full parameter sweep. The returned `sol.original` records which
path succeeded under `:path` (`:trivial` or `:sweep_fallback`).

This is MTK's default `initializealg.nlsolve` when a system contains
`homotopy(...)` nodes. Users override by passing `OverrideInit(nlsolve = ...)`
explicitly.
"""
struct TrivialThenSweep{T, S} <: SciMLBase.AbstractNonlinearAlgorithm
    trivial::T
    sweep::S
end

TrivialThenSweep(; trivial = TrivialHomotopy(), sweep) =
    TrivialThenSweep(trivial, sweep)

function CommonSolve.solve(prob::NonlinearProblem, alg::TrivialThenSweep; kwargs...)
    trivial_sol = CommonSolve.solve(prob, alg.trivial; kwargs...)
    if successful_retcode(trivial_sol)
        return _annotate_path(trivial_sol, :trivial)
    end
    sweep_sol = CommonSolve.solve(prob, alg.sweep; kwargs...)
    return _annotate_path(sweep_sol, :sweep_fallback)
end

# Attach a path marker to sol via the `original` field. The inner solver's
# previous `original` (if any) is preserved under `:inner` so callers that need
# the underlying object can still reach it.
function _annotate_path(sol, path::Symbol)
    inner_original = hasproperty(sol, :original) ? getproperty(sol, :original) : nothing
    return @set sol.original = (; path = path, inner = inner_original)
end
