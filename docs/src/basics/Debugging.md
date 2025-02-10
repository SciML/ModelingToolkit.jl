# Debugging

Every (mortal) modeler writes models that contain mistakes or are susceptible to numerical errors in their hunt for the perfect model.
Debugging such errors is part of the modeling process, and ModelingToolkit includes some functionality that helps with this.

For example, consider an ODE model with "dangerous" functions (here `√`):

```@example debug
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D

@variables u1(t) u2(t) u3(t)
eqs = [D(u1) ~ -√(u1), D(u2) ~ -√(u2), D(u3) ~ -√(u3)]
defaults = [u1 => 1.0, u2 => 2.0, u3 => 3.0]
@named sys = ODESystem(eqs, t; defaults)
sys = structural_simplify(sys)
```

This problem causes the ODE solver to crash:

```@repl debug
prob = ODEProblem(sys, [], (0.0, 10.0), []);
sol = solve(prob, Tsit5());
```

This suggests *that* something went wrong, but not exactly *what* went wrong and *where* it did.
In such situations, the `debug_system` function is helpful:

```@repl debug
dsys = debug_system(sys; functions = [sqrt]);
dprob = ODEProblem(dsys, [], (0.0, 10.0), []);
dsol = solve(dprob, Tsit5());
```

Now we see that it crashed because `u1` decreased so much that it became negative and outside the domain of the `√` function.
We could have figured that out ourselves, but it is not always so obvious for more complex models.

Suppose we also want to validate that `u1 + u2 >= 2.0`. We can do this via the assertions functionality.

```@example debug
@mtkbuild sys = ODESystem(eqs, t; defaults, assertions = [(u1 + u2 >= 2.0) => "Oh no!"])
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
dsol = solve(dprob, Tsit5());
```

Note the logs containing the failed assertion and corresponding message. To temporarily disable
logging in a system returned from `debug_system`, use `ModelingToolkit.ASSERTION_LOG_VARIABLE`.

```@repl debug
dprob[ModelingToolkit.ASSERTION_LOG_VARIABLE] = false;
solve(drob, Tsit5());
```

```@docs
debug_system
```
