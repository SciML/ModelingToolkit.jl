# [Multi- and Nonlocal- Continuation](@id attractors)

In the tutorial on [Bifurcation Diagrams](@ref bifurcation_diagrams) we saw how one can create them by integrating ModelingToolkit.jl with BifurcationKit.jl.
This approach is also often called _continuation_ in the broader literature,
because in essence we are "continuing" the location of individual un/stable fixed points or limit cycles in a dynamical system across a parameter axis.

Recently, an alternative continuation framework was proposed that takes a fundamentally different approach to continuation that is particularly suitable for complex systems. This framework is implemented in [Attractors.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/attractors/stable/) as part of the DynamicalSystems.jl software library.
This new continuation is called _global_ continuation, while the one of BifurcationKit.jl is called _local_ continuation.

Instead of continuing an individual fixed point or limit cycle, the global continuation finds all attractors of the dynamical system and continues all of them, in parallel, in a single continuation. It distinguishes and labels automatically the different attractors.
Hence "multi-" for multiple attractors.
Another key difference is that instead of estimating the local (or linear, or Jacobian) stability of the attractors, it estimates various measures of _nonlocal_ stability (e.g, related with the size of the basins of attraction, or the size of a perturbation that would make the dynamical system state converge to an alternative attractor).
Hence the "nonlocal-" component.
More differences and pros & cons are discussed in the documentation of Attractors.jl.

!!! note "Attractors and basins"
    
    This tutorial assumes that you have some familiarity with dynamical systems,
    specifically what are attractors and basins of attraction. If you don't have
    this yet, we recommend Chapter 1 of the textbook
    [Nonlinear Dynamics](https://link.springer.com/book/10.1007/978-3-030-91032-7).

## Creating the `DynamicalSystem` via MTK

Let's showcase this framework by modelling a chaotic bistable dynamical system that we define via ModelingToolkit.jl, which will the be casted into a `DynamicalSystem` type for the DynamicalSystems.jl library. The equations of our system are

```@example Attractors
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@variables x(t)=-4.0 y(t)=5.0 z(t)=0.0
@parameters a=5.0 b=0.1

eqs = [
    D(x) ~ y - x,
    D(y) ~ -x * z + b * abs(z),
    D(z) ~ x * y - a
]
```

Because our dynamical system is super simple, we will directly make an `System` and cast it in an `ODEProblem` as in the [`Systems` tutorial](@ref programmatically). Since all state variables and parameters have a default value we can immediately write

```@example Attractors
@named modlorenz = System(eqs, t)
ssys = mtkcompile(modlorenz)
# The timespan given to the problem is irrelevant for DynamicalSystems.jl
prob = ODEProblem(ssys, [], (0.0, 1.0))
```

This `prob` can be turned to a dynamical system as simply as

```@example Attractors
using Attractors # or `DynamicalSystems`
ds = CoupledODEs(prob)
```

DynamicalSystems.jl integrates fully with ModelingToolkit.jl out of the box and understands whether a given problem has been created via ModelingToolkit.jl. For example you can use the symbolic variables, or their `Symbol` representations, to access a system state or parameter

```@example Attractors
observe_state(ds, x)
```

```@example Attractors
current_parameter(ds, :a) # or `a` directly
```

## Finding all attractors in the state space

Attractors.jl provides an extensive interface for finding all (within a state space region and numerical accuracy) attractors of a dynamical system.
This interface is structured around the type `AttractorMapper` and is discussed in the Attractors.jl documentation in detail. Here we will briefly mention one of the possible approaches, the recurrences-based algorithm. It finds attractors by finding locations in the state space where the trajectory returns again and again.

To use this technique, we first need to create a tessellation of the state space

```@example Attractors
grid = (
    range(-15.0, 15.0; length = 150), # x
    range(-20.0, 20.0; length = 150), # y
    range(-20.0, 20.0; length = 150) # z
)
```

which we then give as input to the `AttractorsViaRecurrences` mapper along with the dynamical system

```@example Attractors
mapper = AttractorsViaRecurrences(ds, grid;
    consecutive_recurrences = 1000,
    consecutive_lost_steps = 100
)
```

to learn about the metaparameters of the algorithm visit the documentation of Attractors.jl.

This `mapper` object is incredibly powerful! It can be used to map initial conditions to attractor they converge to, while ensuring that initial conditions that converge to the same attractor are given the same label.
For example, if we use the `mapper` as a function and give it an initial condition we get

```@example Attractors
mapper([-4.0, 5, 0])
```

```@example Attractors
mapper([4.0, 2, 0])
```

```@example Attractors
mapper([1.0, 3, 2])
```

The numbers returned are simply the unique identifiers of the attractors the initial conditions converged to.

DynamicalSystems.jl library is the only dynamical systems software (in any language) that provides such an infrastructure for mapping initial conditions of any arbitrary dynamical system to its unique attractors. And this is only the tip of this iceberg! The rest of the functionality of Attractors.jl is all full of brand new cutting edge progress in dynamical systems research.

The found attractors are stored in the mapper internally, to obtain them we
use the function

```@example Attractors
attractors = extract_attractors(mapper)
```

This is a dictionary that maps attractor IDs to the attractor sets themselves.
`StateSpaceSet` is a wrapper of a vector of points and behaves exactly like a vector of points. We can plot them easily like

```@example Attractors
using CairoMakie
fig = Figure()
ax = Axis(fig[1, 1])
colors = ["#7143E0", "#191E44"]
for (id, A) in attractors
    scatter!(ax, A[:, [1, 3]]; color = colors[id])
end
fig
```

## Basins of attraction

Estimating the basins of attraction of these attractors is a matter of a couple lines of code.
First we define the state space are to estimate the basins for.
Here we can re-use the `grid` we defined above. Then we only have to call

```julia
basins = basins_of_attraction(mapper, grid)
```

We won't run this in this tutorial because it is a length computation (150×150×150).
We will however estimate a slice of the 3D basins of attraction.
DynamicalSystems.jl allows for a rather straightforward setting of initial conditions:

```@example Attractors
ics = [Dict(:x => x, :y => 0, :z => z) for x in grid[1] for z in grid[3]]
```

now we can estimate the basins of attraction on a slice on the x-z grid

```@example Attractors
fs, labels = basins_fractions(mapper, ics)
labels = reshape(labels, (length(grid[1]), length(grid[3])))
```

and visualize them

```@example Attractors
heatmap(grid[1], grid[3], labels; colormap = colors)
```

## Global continuation

We've already outlined the principles of the global continuation, so let's just do it here!
We first have to define a global continuation algorithm, which for this tutorial,
it is just a wrapper of the existing `mapper`

```@example Attractors
ascm = AttractorSeedContinueMatch(mapper);
```

we need two more ingredients to perform the global continuation.
One is a sampler of initial conditions in the state space.
Here we'll uniformly sample initial conditions within this grid we have already defined

```@example Attractors
sampler, = statespace_sampler(grid);
```

the last ingredient is what parameter(s) to perform the continuation over.
In contrast to local continuation, where we can only specify a parameter range, in global continuation one can specify an exact parameter curve to continue over.
This curve can span any-dimensional parameter space, in contrast to the 1D or 2D parameter spaces supported in local continuation.
Here we will use the curve

```@example Attractors
params(θ) = [:a => 5 + 0.5cos(θ), :b => 0.1 + 0.01sin(θ)]
θs = range(0, 2π; length = 101)
pcurve = params.(θs)
```

which makes an ellipsis over the parameter space.

We put these three ingredients together to call the global continuation

```@example Attractors
fractions_cont, attractors_cont = global_continuation(ascm, pcurve, sampler);
```

The output of the continuation is how the attractors and their basins fractions change over this parameter curve. We can visualize this directly using a convenience function

```@example Attractors
fig = plot_basins_attractors_curves(
    fractions_cont, attractors_cont, A -> minimum(A[:, 1]), θs
)
```

The top panel shows the relative basins of attractions of the attractors and the bottom panel shows their minimum x-position. The colors correspond to unique attractors. Perhaps making a video is easier to understand:

```@example Attractors
animate_attractors_continuation(
    ds, attractors_cont, fractions_cont, pcurve;
    savename = "curvecont.mp4"
);
```

```@raw html
<video width="auto" controls loop>
<source src="../curvecont.mp4" type="video/mp4">
</video>
```

To learn more about this global continuation and its various options, and more details about how it compares with local continuation, visit the documentation of [Attractors.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/attractors/stable/).
