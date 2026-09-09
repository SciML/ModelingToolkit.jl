# Checkpoint: array equations through mtkcompile / ODEProblem

**Status:** partial / work-in-progress checkpoint opened because agent usage was about to run out.
**Date:** 2026-09-09
**Related PRs:**
- ModelingToolkit (this PR / `array-aware-tearing`): preserve through `mtkcompile` via array equation groups
- ModelingToolkit [#5101](https://github.com/SciML/ModelingToolkit.jl/pull/5101): `complete` → `ODEProblem` for array DEs (codegen only; no mtkcompile)
- ModelingToolkit [#5100](https://github.com/SciML/ModelingToolkit.jl/pull/5100): **closed** — `scalarize_arrays=false` / skip tearing was rejected as a non-useful intermediate
- StateSelection / ModelingToolkitTearing: companion PR `mtkt-array-equation-groups` (ArrayEquationGroup tracking)
- MethodOfLines [#688](https://github.com/SciML/MethodOfLines.jl/pull/688): still based on the closed `scalarize_arrays=false` story; needs rewrite after MTK lands

## Goal (Chris)

1. Array differential equations should build an `ODEProblem` (and DAEProblem) without requiring scalarization of the generated code (O(1) in array length).
2. The **real** hard problem: index reduction and tearing must **keep array equations intact**, not turn those algorithms off.

MethodOfLines today gets O(1) only by skipping `mtkcompile` and using `DAEProblem` from `complete`.

## What is done

### A. `complete` → `ODEProblem` (PR #5101)

On a completed system (no `mtkcompile`), array DEs of the form `D(u[slice]) ~ f` or residual `D(u[slice]) - f ~ 0` can build an `ODEProblem`. Codegen packs into `du` via the same ArrayMaker / contiguous-view ideas as the DAE residual path. This unstacks from #5100 and does **not** require a new keyword.

**Verify:** tests in `lib/ModelingToolkitBase/test/array_equation_ode.jl` (on the #5101 branch).

### B. Array equation groups through structural simplification (this checkpoint + StateSelection)

Design choice (not "turn tearing off"):

- Array eqs are still **scalarized into bipartite graph rows** so matching, Pantelides, dummy derivatives, tearing, and alias elimination see exact per-element incidence.
- Rows remember the parent array equation (`ArrayEquationGroup`, `row_group`, `row_elem`).
- Passes that rewrite a row in a way the array equation cannot represent mark the group **dirty** (differentiation of the eq, removal, dummy-derivative substitution, solving for another variable, inline linear SCCs, clock partition splits, …).
- Intact groups are **excluded from integer-linear Gaussian elimination** as pivots/reducees (`linear_subsys_adjmat!` / `is_intact_array_group_row`) so alias elimination cannot break the group.
- With `preserve_array_equations = true` on `DefaultReassembleAlgorithm` / `mtkcompile`, intact groups are **reassembled** into a single array equation over scalar unknowns; dirty groups stay scalarized.
- Default `mtkcompile` output is unchanged (`preserve_array_equations = false`).

MTK side of this checkpoint also has:
- docs in `docs/src/internals/mtkcompile.md`
- tests in `test/structural_transformation/array_equations.jl`
- codegen helpers for O(1) array literals / views when emitting preserved eqs
- alias-elimination awareness of groups

**MTKTearing unit tests:** local run of `ModelingToolkitTearing/test/runtests.jl` was all Pass.
**Broader MTK regression sample** (tearing, array_equations, odesystem, init, codegen, SII): filtered Pass after fixing test env deps.

## What is NOT done / required next

### Must do next (hard path)

1. **Atomic array incidence (if Chris requires it):** today's design still expands to scalar rows for the graph. If "keep array equations intact" means the bipartite graph itself must treat an array DE as one node (true O(1) structural simplify, not just O(1) codegen), that is still unimplemented. Dirty-tracking + reassemble is an intermediate that keeps algorithms correct on scalars while restoring array form at the end when safe.
2. **Wire `preserve_array_equations` end-to-end on released MTK:** ensure `mtkcompile(...; preserve_array_equations=true)` forwards into `DefaultReassembleAlgorithm` on the versions used in CI; bump ModelingToolkitTearing compat once StateSelection PR is merged/registered.
3. **Default-on policy:** decide when preserve can be the default for eligible eqs without a kwarg (rollout safety vs Chris's "prefer default-correct").
4. **Dirty-marking audit:** prove we dirty on every pass that breaks representability, and do **not** dirty on benign substitutions (observed aliases). Mix of array DE + nonlinear algebraics + high-index systems needs more tests.
5. **High-index / must-scalarize cases:** Pantelides that differentiates an array eq currently dirties → scalar emit. Confirm that is correct; if an entire array eq can be differentiated as a unit and stay an array, implement that instead of dirtying.
6. **ODE/DAE codegen after preserve-`mtkcompile`:** #5101 covers `complete` only. After reassemble, `ODEProblem`/`DAEProblem` must accept the preserved array eqs (may already work if #5101 checks are shape-based rather than `arrays_scalarized` metadata). Add integration tests: `mtkcompile(sys; preserve_array_equations=true)` → `ODEProblem` / `DAEProblem`, Expr size independent of `n`.
7. **Initialization:** #5093 (`vec(::Bool)` for array `uˍt` guess) and related init path; MOL still uses `build_initializeprob=false` workarounds.
8. **#5097:** `expand_array_derivatives!` still tends to unroll into `array_literal` of scalar `du` reads on some paths — prefer contiguous `view(du, …)` for O(1) AST.
9. **MethodOfLines:** rewrite #688 off `scalarize_arrays=false`; use `complete`→ODE (#5101) and/or preserve-`mtkcompile` once available; assert O(1) treesize for heat-style ODE path like the DAE suite.
10. **Performance:** structural simplify is still O(N) in graph size even when codegen is O(1). True O(1) compile needs atomic array nodes or compressed incidence.

### Explicitly rejected

- `mtkcompile(; scalarize_arrays=false)` / `__mtkcompile_no_tearing` as the solution (#5100 closed).

## How to verify locally

```julia
# 1) complete → ODE (needs #5101)
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
n = 96
@variables u(t)[1:n]
sys = complete(System([D(u) ~ -u], t; name=:m))
prob = ODEProblem(sys, [u => ones(n)], (0.0, 1.0))
solve(prob, Tsit5())

# 2) preserve through mtkcompile (needs this PR + MTKTearing from StateSelection branch)
# sys = mtkcompile(heat; preserve_array_equations=true)
# @assert length(equations(sys)) == 1   # for pure interior array DE + observed BCs
```

## Suggested PR stack when usage resumes

1. Land / fix CI on #5101 (`complete`→ODE).
2. Land StateSelection `mtkt-array-equation-groups`.
3. Land MTK `array-aware-tearing` against a MTKTearing version that includes (2); add end-to-end O(1) tests.
4. Rewrite MethodOfLines #688.
5. (Larger follow-up) Atomic array nodes in the bipartite graph for true O(1) structural simplify.
