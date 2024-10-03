# [Automatic Transformation of Nth Order ODEs to 1st Order ODEs](@id higher_order_ODE)

ModelingToolkit has a system for transformations of mathematical
systems. These transformations allow for symbolically changing
the representation of the model to problems that are easier to
numerically solve. One simple to demonstrate transformation is the
`structural_simplify`, which does a lot of tricks, one being the
transformation that turns an Nth order ODE into N
coupled 1st order ODEs.

To see this, let's define a second order riff on the Lorenz equations.
We utilize the derivative operator twice here to define the second order:

```@example orderlowering
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters σ ρ β
@variables x(t) y(t) z(t)

eqs = [D(D(x)) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z]

@named sys = ODESystem(eqs, t)
```

Note that we could've used an alternative syntax for 2nd order, i.e.
`D = Differential(t)^2` and then `D(x)` would be the second derivative,
and this syntax extends to `N`-th order. Also, we can use `*` or `∘` to compose
`Differential`s, like `Differential(t) * Differential(x)`.

Now let's transform this into the `ODESystem` of first order components.
We do this by calling `structural_simplify`:

```@example orderlowering
sys = structural_simplify(sys)
```

Now we can directly numerically solve the lowered system. Note that,
following the original problem, the solution requires knowing the
initial condition for `x'`, and thus we include that in our input
specification:

```@example orderlowering
u0 = [D(x) => 2.0,
    x => 1.0,
    y => 0.0,
    z => 0.0]

p = [σ => 28.0,
    ρ => 10.0,
    β => 8 / 3]

tspan = (0.0, 100.0)
prob = ODEProblem(sys, u0, tspan, p, jac = true)
sol = solve(prob, Tsit5())
using Plots;
plot(sol, idxs = (x, y));
```
