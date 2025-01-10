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

```@docs
debug_system
```
