# [Homotopy](@id homotopy)

ModelingToolkit implements the Modelica `homotopy(actual, simplified)` operator
([Modelica Specification 3.7.4.2](https://specification.modelica.org/master/operators-and-expressions.html#homotopy))
as a way to robustly solve nonlinear systems that are hard to solve from a cold
start. It is most commonly reached for during initialization, but it is a general
nonlinear-solving construct: any nonlinear system whose equations carry a
`homotopy` annotation can be solved by continuation.

## What Is Homotopy?

The `simplified` equations are a set of equations for which the nonlinear system
is easier to get a convergent (Newton) iteration for, and the `actual` equations
are the more complex equations you actually want to solve. A homotopy solver
starts by solving the nonlinear system with the `simplified` equations, and uses
that solution as the starting point for solving the `actual` problem, deforming
continuously from one to the other.

There are many ways this can be used. For example, if you have equations with
multiple solutions — like a quadratic equation with a positive and a negative
root — you can simplify it down to an approximating linear problem that has a
single (say, positive) solution, and then deform it to the actual equation to
stabilise the process of converging to that positive solution.

The operator encodes both expressions in a single annotation:

```julia
homotopy(actual, simplified)
```

Concretely, the continuation introduces a scalar parameter ``\lambda`` and solves

```math
(1 - \lambda)\,\text{simplified} + \lambda\,\text{actual}
```

sweeping ``\lambda`` from `0` (the easy `simplified` system) to `1` (the
`actual` system), warm-starting each step from the previous solution.

## Runtime Semantics

The operator stays an opaque symbolic function through `System` construction,
`mtkcompile`, and runtime code generation — no continuation parameter is added to
the system. Wherever the operator is evaluated numerically outside a continuation
solve, the generated code calls the numeric fallback

```julia
homotopy(actual::Real, simplified::Real) = actual
```

so the operator evaluates to `actual`, as the Modelica specification prescribes.
Note the honest cost: the `simplified` argument expression is still evaluated and
its value discarded (arguments are evaluated before the call). This small
overhead is borne only by systems that use `homotopy` — systems without the
operator go through a byte-identical pipeline and are completely unaffected.

Symbolic differentiation works through the operator: nodewise derivative rules
keep symbolic jacobians, `tgrad`, and index reduction consistent. At runtime,
differentiated equations reproduce `actual`'s derivative; along the continuation
they follow the derivative of the blended expression above.

## Building a `HomotopyProblem`

A system whose equations contain `homotopy` nodes is built into a
[`SciMLBase.HomotopyProblem`](@ref) whose residual is the blended expression
above, compiled as `f(u, p, λ)`. `λ` is an explicit trailing argument — it is
never added to the system's parameters, and your parameter object `p` passes
through untouched. All `homotopy` calls in a system share the single `λ`, per the
Modelica spec's recommendation of (conceptually) one homotopy iteration over the
whole model; this includes nested `homotopy` calls.

There are two ways to construct it:

```julia
# Explicit: always returns a HomotopyProblem (errors if `sys` has no `homotopy`).
prob = HomotopyProblem(sys, op)

# Automatic: returns a HomotopyProblem when `sys` contains `homotopy` nodes, and
# a plain NonlinearProblem otherwise.
prob = AbstractNonlinearProblem(sys, op)
```

The `HomotopyProblem`'s `λspan` defaults to `(0.0, 1.0)`. It is solved by
`NonlinearSolveBase.HomotopySweep`, a natural-parameter continuation solver that
sweeps `λ` from `0` to `1`, solving a standard nonlinear problem at each step and
warm-starting from the previous step's solution. A `HomotopyProblem` with no
algorithm (`solve(prob)`) defaults to this solver.

## Example: Out-of-Basin Rescue

The equation `0 = atan(y - 3)` has a root at `y = 3`, but a Newton solver
starting from `y = 12` diverges because `atan` saturates. Using `homotopy` with
`simplified = y` (whose root is `y = 0`) lets the continuation walk from the easy
root to the true one:

```julia
using ModelingToolkit, NonlinearSolve

@variables y
@mtkcompile sys = System([0 ~ homotopy(atan(y - 3), y)])
prob = HomotopyProblem(sys, [y => 12.0])
sol = solve(prob, HomotopySweep())
sol[y] # ≈ 3.0 — the continuation rescued the out-of-basin guess
```

The operating point (`[y => 12.0]`) provides the starting point of the
continuation; the sweep deforms the equations so the solver reaches `y ≈ 3` at
`λ = 1`.

## Broadcasting Over Arrays

`homotopy` is a scalar operator. For array equations, broadcast it elementwise:

```julia
eqs = 0 .~ homotopy.(actual_array, simplified_array)
```

This creates one `homotopy` node per element; the continuation lowering rewrites
each node independently, and all of them share the single continuation parameter
`λ`.

## Customizing the Continuation Solver

To tune the sweep, pass your own continuation algorithm to `solve`:

```julia
sol = solve(prob, HomotopySweep(nsteps = 30))
```

`HomotopySweep` accepts the keyword arguments `inner` (the nonlinear algorithm
used at each step), `nsteps`, `adaptive`, `initial_step_factor`, and `min_dλ`;
see the `NonlinearSolveBase.HomotopySweep` docstring for their meanings and
defaults.

## Limitations

  - **`expression = Val{true}` is not yet supported** for the homotopy
    constructor; build the problem directly (the default
    `expression = Val{false}`). This can be added in a future PR.
  - **The jacobian/sparsity of the standard build are dropped.** They encode the
    `λ = 1` (opaque-`actual`) system and would be wrong mid-sweep; continuation
    steps solve with a freshly differentiated residual. Per-problem analytic
    jacobians for the swept residual are future work.
  - **Scalar `Real` expressions only** for a single `homotopy` node, matching
    Modelica's restriction; use broadcasting (above) for arrays.
  - **Only equations and observed equations are rewritten.** A `homotopy` call
    inside a parameter binding or default value is left as-is and evaluates as
    `actual` at all `λ`.

## API Reference

```@docs
ModelingToolkit.homotopy
```

## See Also

  - [Modelica Specification 3.7.4.2](https://specification.modelica.org/master/operators-and-expressions.html#homotopy)
    — the upstream specification this operator implements.
