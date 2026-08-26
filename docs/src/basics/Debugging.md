# Debugging

Every (mortal) modeler writes models that contain mistakes or are susceptible to numerical errors in their hunt for the perfect model.
Debugging such errors is part of the modeling process, and ModelingToolkit includes some functionality that helps with this.

For example, consider an ODE model with "dangerous" functions (here `√`):

```@example debug
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D

@variables u1(t) u2(t) u3(t)
eqs = [D(u1) ~ -√(u1), D(u2) ~ -√(u2), D(u3) ~ -√(u3)]
initial_conditions = [u1 => 1.0, u2 => 2.0, u3 => 3.0]
@named sys = System(eqs, t; initial_conditions)
sys = mtkcompile(sys)
```

This problem causes the ODE solver to crash:

```@repl debug
prob = ODEProblem(sys, [], (0.0, 10.0));
sol = solve(prob, Tsit5());
```

This suggests *that* something went wrong, but not exactly *what* went wrong and *where* it did.
In such situations, the `debug_system` function is helpful:

```@repl debug
dsys = debug_system(sys; functions = [sqrt]);
dprob = ODEProblem(dsys, [], (0.0, 10.0));
dsol = solve(dprob, Tsit5());
```

Now we see that it crashed because `u1` decreased so much that it became negative and outside the domain of the `√` function.
We could have figured that out ourselves, but it is not always so obvious for more complex models.

Suppose we also want to validate that `u1 + u2 >= 2.0`. We can do this via the assertions functionality.

```@example debug
@mtkcompile sys = System(eqs, t; initial_conditions, assertions = [(u1 + u2 >= 2.0) => "Oh no!"])
```

The assertions must be an iterable of pairs, where the first element is the symbolic condition and
the second is a message to be logged when the condition fails. All assertions are added to the
generated code and will cause the solver to reject steps that fail the assertions. For systems such
as the above where the assertion is guaranteed to eventually fail, the solver will likely exit
with a `dtmin` failure..

```@example debug
prob = ODEProblem(sys, [], (0.0, 10.0))
sol = solve(prob, Tsit5())
```

We can use `debug_system` to log the failing assertions in each call to the RHS function.

```@repl debug
dsys = debug_system(sys; functions = []);
dprob = ODEProblem(dsys, [], (0.0, 10.0));
dsol = solve(dprob, Tsit5(); dtmin = 0.1); # high dtmin only to show less clutter on this page
```

Note the logs containing the failed assertion and corresponding message. To temporarily disable
logging in a system returned from `debug_system`, use `ModelingToolkit.ASSERTION_LOG_VARIABLE`.

```@repl debug
dprob[ModelingToolkit.ASSERTION_LOG_VARIABLE] = false;
solve(dprob, Tsit5());
```

```@docs; canonical = false
debug_system
```

## Controlling verbosity

Diagnostic output from ModelingToolkit — the `mtkcompile` pipeline, problem construction
and initialization, and analysis utilities — is controlled through
[SciMLLogging.jl](https://github.com/SciML/SciMLLogging.jl) via the `verbose` keyword
argument, which accepts an [`MTKVerbosity`](@ref) specifier, a `SciMLLogging` preset, or
a `Bool`:

```julia
using SciMLLogging: SciMLLogging, Silent, InfoLevel

mtkcompile(sys)                                    # default (Standard) verbosity
mtkcompile(sys; verbose = false)                   # silence everything
mtkcompile(sys; verbose = SciMLLogging.Detailed()) # more diagnostics

# fine-grained control per toggle or group
mtkcompile(sys; verbose = MTKVerbosity(state_priority_tie = Silent))
mtkcompile(sys; verbose = MTKVerbosity(compilation = InfoLevel))
```

The available toggles are documented on [`MTKVerbosity`](@ref):

| Toggle                            | Default     | Controls                                                               |
|:--------------------------------- |:----------- |:---------------------------------------------------------------------- |
| `state_priority_tie`              | `WarnLevel` | Tied `state_priority` in an alias group during alias elimination        |
| `underconstrained_variables`      | `Silent`    | Report of underconstrained variables found during alias elimination     |
| `if_lifting_condition_grammar`    | `WarnLevel` | `IfLifting` skipping a condition outside its supported grammar          |
| `observed_equation_cycle`         | `Silent`    | Printing the smallest cycle when sorting observed equations fails       |
| `singular_initialization`         | `WarnLevel` | Structurally singular initialization system                             |
| `overdetermined_initialization`   | `WarnLevel` | Overdetermined initialization system (least-squares fallback)           |
| `underdetermined_initialization`  | `WarnLevel` | Underdetermined initialization system (least-squares fallback)          |
| `scc_initialization_unavailable`  | `WarnLevel` | `SCCNonlinearProblem` initialization requires `split = true`            |
| `cyclic_dependency`               | `Silent`    | Cycles among initial conditions of unknowns/parameters                  |
| `overdetermined_constraints`      | `WarnLevel` | Overdetermined `BVProblem`/dynamic-optimization constraints             |
| `missing_scc_schedule`            | `WarnLevel` | Simplified system unexpectedly missing a tearing schedule               |
| `dynamic_opt_time_grid`           | `WarnLevel` | Ignored `dt`/`steps` argument in dynamic optimization                   |
| `empty_operating_point`           | `WarnLevel` | Empty operating point passed to `linearization_function`                |
| `initialization_analysis`         | `InfoLevel` | The report of `analyze_initialization_jacobian`                         |
| `no_unbound_inputs`               | `WarnLevel` | `generate_control_function` found no unbound inputs                     |
| `analysis_point_causality`        | `WarnLevel` | Reversed causality in an analysis-point `connect`                       |

### Problem constructors

`SciMLBase.*Problem` constructors share the `verbose` keyword with the solver, since
problem keyword arguments are forwarded to `solve`. The value is routed by type: an
`MTKVerbosity` is consumed by ModelingToolkit and *not* placed in `prob.kwargs`; a
preset or `Bool` applies to ModelingToolkit *and* is forwarded to the solver; a solver
verbosity specifier (e.g. `DEVerbosity`) is forwarded untouched. The
`initialization_verbosity` sub-specifier of `MTKVerbosity` controls the verbosity of the
initialization problem's solve (default `Minimal()`), and accepts a preset or a
`NonlinearVerbosity`.

The boolean keywords `warn_initialize_determined`, `warn_cyclic_dependency`, and
`warn_empty_op` are deprecated in favor of `verbose`; when explicitly passed they
override the corresponding toggles.

```@docs; canonical = false
MTKVerbosity
```

## Error message guidance

Errors raised while building a problem describe what went wrong and then, where it is
useful, how to fix it in terms of the ModelingToolkit API — which keyword argument to pass,
or where in the system to supply a value. Front ends which present a modelling language of
their own can leave that part out, since advice to call `ODEProblem` or to pass
`initialization_eqs` would send their users looking for something their language does not
have. The description of what went wrong is unaffected.

```julia
ModelingToolkit.show_api_guidance!(false)
```

!!! warning "Experimental"

    `show_api_guidance!` is experimental and unsupported. It may change or be removed in
    any release, without a breaking version bump. Verbosity across the SciML ecosystem is
    moving to [SciMLLogging.jl](https://github.com/SciML/SciMLLogging.jl), which
    ModelingToolkit has begun adopting via [`MTKVerbosity`](@ref); this setting is
    expected to be absorbed into that mechanism. That will be the supported way to
    control this.
