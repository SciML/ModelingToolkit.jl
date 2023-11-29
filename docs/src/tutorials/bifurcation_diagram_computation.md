# [Bifurcation Diagrams](@id bifurcation_diagrams)

Bifurcation diagrams describes how, for a dynamic system, the quantity and quality of its steady states changes with a parameter's value. These can be computed through the [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl) package. ModelingToolkit provides a simple interface for creating BifurcationKit compatible `BifurcationProblem`s from `NonlinearSystem`s and `ODESystem`s. All the features provided by BifurcationKit can then be applied to these systems. This tutorial provides a brief introduction for these features, with BifurcationKit.jl providing [a more extensive documentation](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/).

### Creating a `BifurcationProblem`

Let us first consider a simple `NonlinearSystem`:

```@example Bif1
using ModelingToolkit
@variables t x(t) y(t)
@parameters μ α
eqs = [0 ~ μ * x - x^3 + α * y,
    0 ~ -y]
@named nsys = NonlinearSystem(eqs, [x, y], [μ, α])
```

we wish to compute a bifurcation diagram for this system as we vary the parameter `μ`. For this, we need to provide the following information:

 1. The system for which we wish to compute the bifurcation diagram (`nsys`).
 2. The parameter which we wish to vary (`μ`).
 3. The parameter set for which we want to compute the bifurcation diagram.
 4. An initial guess of the state of the system for which there is a steady state at our provided parameter value.
 5. The variable which value we wish to plot in the bifurcation diagram (this argument is optional, if not provided, BifurcationKit default plot functions are used).

We declare this additional information:

```@example Bif1
bif_par = μ
p_start = [μ => -1.0, α => 1.0]
u0_guess = [x => 1.0, y => 1.0]
plot_var = x;
```

For the initial state guess (`u0_guess`), typically any value can be provided, however, read BifurcatioKit's documentation for more details.

We can now create our `BifurcationProblem`, which can be provided as input to BifurcationKit's various functions.

```@example Bif1
using BifurcationKit
bprob = BifurcationProblem(nsys,
    u0_guess,
    p_start,
    bif_par;
    plot_var = plot_var,
    jac = false)
```

Here, the `jac` argument (by default set to `true`) sets whenever to provide BifurcationKit with a Jacobian or not.

### Computing a bifurcation diagram

Let us consider the `BifurcationProblem` from the last section. If we wish to compute the corresponding bifurcation diagram we must first declare various settings used by BifurcationKit to compute the diagram. These are stored in a `ContinuationPar` structure (which also contain a `NewtonPar` structure).

```@example Bif1
p_span = (-4.0, 6.0)
opt_newton = NewtonPar(tol = 1e-9, max_iterations = 20)
opts_br = ContinuationPar(dsmin = 0.001, dsmax = 0.05, ds = 0.01,
    max_steps = 100, nev = 2, newton_options = opt_newton,
    p_min = p_span[1], p_max = p_span[2],
    detect_bifurcation = 3, n_inversion = 4, tol_bisection_eigenvalue = 1e-8,
    dsmin_bisection = 1e-9);
```

Here, `p_span` sets the interval over which we wish to compute the diagram.

Next, we can use this as input to our bifurcation diagram, and then plot it.

```@example Bif1
bf = bifurcationdiagram(bprob, PALC(), 2, (args...) -> opts_br; bothside = true)
```

Here, the value `2` sets how sub-branches of the diagram that BifurcationKit should compute. Generally, for bifurcation diagrams, it is recommended to use the `bothside=true` argument.

```@example Bif1
using Plots
plot(bf;
    putspecialptlegend = false,
    markersize = 2,
    plotfold = false,
    xguide = "μ",
    yguide = "x")
```

Here, the system exhibits a pitchfork bifurcation at *μ=0.0*.

### Using `ODESystem` inputs

It is also possible to use `ODESystem`s (rather than `NonlinearSystem`s) as input to `BifurcationProblem`. Here follows a brief such example.

```@example Bif2
using BifurcationKit, ModelingToolkit, Plots

@variables t x(t) y(t)
@parameters μ
D = Differential(t)
eqs = [D(x) ~ μ * x - y - x * (x^2 + y^2),
    D(y) ~ x + μ * y - y * (x^2 + y^2)]
@named osys = ODESystem(eqs, t)

bif_par = μ
plot_var = x
p_start = [μ => 1.0]
u0_guess = [x => 0.0, y => 0.0]

bprob = BifurcationProblem(osys,
    u0_guess,
    p_start,
    bif_par;
    plot_var = plot_var,
    jac = false)

p_span = (-3.0, 3.0)
opt_newton = NewtonPar(tol = 1e-9, max_iterations = 20)
opts_br = ContinuationPar(dsmin = 0.001, dsmax = 0.05, ds = 0.01,
    max_steps = 100, nev = 2, newton_options = opt_newton,
    p_max = p_span[2], p_min = p_span[1],
    detect_bifurcation = 3, n_inversion = 4, tol_bisection_eigenvalue = 1e-8,
    dsmin_bisection = 1e-9)

bf = bifurcationdiagram(bprob, PALC(), 2, (args...) -> opts_br; bothside = true)
using Plots
plot(bf;
    putspecialptlegend = false,
    markersize = 2,
    plotfold = false,
    xguide = "μ",
    yguide = "x")
```

Here, the value of `x` in the steady state does not change, however, at `μ=0` a Hopf bifurcation occur and the steady state turn unstable.
