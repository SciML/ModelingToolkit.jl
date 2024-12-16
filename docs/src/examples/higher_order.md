# Automatic Transformation of Nth Order ODEs to 1st Order ODEs

ModelingToolkit has a system for transformations of mathematical
systems. These transformations allow for symbolically changing
the representation of the model to problems that are easier to
numerically solve. One simple to demonstrate transformation, is
`structural_simplify`, which does a lot of tricks, one being the
transformation that turns an Nth order ODE into N
coupled 1st order ODEs.

To see this, let's define a second order riff on the Lorenz equations.
We utilize the derivative operator twice here to define the second order:

```@example orderlowering
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D

@mtkmodel SECOND_ORDER begin
    @parameters begin
        σ = 28.0
        ρ = 10.0
        β = 8 / 3
    end
    @variables begin
        x(t) = 1.0
        y(t) = 0.0
        z(t) = 0.0
    end
    @equations begin
        D(D(x)) ~ σ * (y - x)
        D(y) ~ x * (ρ - z) - y
        D(z) ~ x * y - β * z
    end
end
@mtkbuild sys = SECOND_ORDER()
```

The second order ODE has been automatically transformed to two first order ODEs.

Note that we could've used an alternative syntax for 2nd order, i.e.
`D = Differential(t)^2` and then `D(x)` would be the second derivative,
and this syntax extends to `N`-th order. Also, we can use `*` or `∘` to compose
`Differential`s, like `Differential(t) * Differential(x)`.

Now let's transform this into the `ODESystem` of first order components.
We do this by calling `structural_simplify`:

Now we can directly numerically solve the lowered system. Note that,
following the original problem, the solution requires knowing the
initial condition for both `x` and `D(x)`.
The former already got assigned a default value in the `@mtkmodel`,
but we still have to provide a value for the latter.

```@example orderlowering
u0 = [D(sys.x) => 2.0]
tspan = (0.0, 100.0)
prob = ODEProblem(sys, u0, tspan, [], jac = true)
sol = solve(prob, Tsit5())
using Plots
plot(sol, idxs = (sys.x, sys.y))
```
