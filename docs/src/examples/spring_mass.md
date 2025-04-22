# Component-Based Modeling of a Spring-Mass System

In this tutorial, we will build a simple component-based model of a spring-mass system. A spring-mass system consists of one or more masses connected by springs. [Hooke's law](https://en.wikipedia.org/wiki/Hooke%27s_law) gives the force exerted by a spring when it is extended or compressed by a given distance. This specifies a differential-equation system where the acceleration of the masses is specified using the forces acting on them.

## Copy-Paste Example

```@example component
using ModelingToolkit, Plots, OrdinaryDiffEq, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D
using Symbolics: scalarize

function Mass(; name, m = 1.0, xy = [0.0, 0.0], u = [0.0, 0.0])
    ps = @parameters m = m
    sts = @variables pos(t)[1:2]=xy v(t)[1:2]=u
    eqs = scalarize(D.(pos) .~ v)
    ODESystem(eqs, t, [pos..., v...], ps; name)
end

function Spring(; name, k = 1e4, l = 1.0)
    ps = @parameters k=k l=l
    @variables x(t), dir(t)[1:2]
    ODESystem(Equation[], t, [x, dir...], ps; name)
end

function connect_spring(spring, a, b)
    [spring.x ~ norm(scalarize(a .- b))
     scalarize(spring.dir .~ scalarize(a .- b))]
end

function spring_force(spring)
    -spring.k .* scalarize(spring.dir) .* (spring.x - spring.l) ./ spring.x
end

m = 1.0
xy = [1.0, -1.0]
k = 1e4
l = 1.0
center = [0.0, 0.0]
g = [0.0, -9.81]
@named mass = Mass(m = m, xy = xy)
@named spring = Spring(k = k, l = l)

eqs = [connect_spring(spring, mass.pos, center)
       scalarize(D.(mass.v) .~ spring_force(spring) / mass.m .+ g)]

@named _model = ODESystem(eqs, t, [spring.x; spring.dir; mass.pos], [])
@named model = compose(_model, mass, spring)
sys = structural_simplify(model)

prob = ODEProblem(sys, [], (0.0, 3.0))
sol = solve(prob, Rosenbrock23())
plot(sol)
```

## Explanation

### Building the components

For each component, we use a Julia function that returns an `ODESystem`. At the top, we define the fundamental properties of a `Mass`: it has a mass `m`, a position `pos` and a velocity `v`. We also define that the velocity is the rate of change of position with respect to time.

```@example component
function Mass(; name, m = 1.0, xy = [0.0, 0.0], u = [0.0, 0.0])
    ps = @parameters m = m
    sts = @variables pos(t)[1:2]=xy v(t)[1:2]=u
    eqs = scalarize(D.(pos) .~ v)
    ODESystem(eqs, t, [pos..., v...], ps; name)
end
```

Note that this is an incompletely specified `ODESystem`. It cannot be simulated on its own, since the equations for the velocity `v[1:2](t)` are unknown. Notice the addition of a `name` keyword. This allows us to generate different masses with different names. A `Mass` can now be constructed as:

```@example component
Mass(name = :mass1)
```

Or using the `@named` helper macro

```@example component
@named mass1 = Mass()
```

Next, we build the spring component. It is characterized by the spring constant `k` and the length `l` of the spring when no force is applied to it. The state of a spring is defined by its current length and direction.

```@example component
function Spring(; name, k = 1e4, l = 1.0)
    ps = @parameters k=k l=l
    @variables x(t), dir(t)[1:2]
    ODESystem(Equation[], t, [x, dir...], ps; name)
end
```

We now define functions that help construct the equations for a mass-spring system. First, the `connect_spring` function connects a `spring` between two positions `a` and `b`. Note that `a` and `b` can be the `pos` of a `Mass`, or just a fixed position such as `[0., 0.]`. In that sense, the length of the spring `x` is given by the length of the vector `dir` joining `a` and `b`.

```@example component
function connect_spring(spring, a, b)
    [spring.x ~ norm(scalarize(a .- b))
     scalarize(spring.dir .~ scalarize(a .- b))]
end
```

Lastly, we define the `spring_force` function that takes a `spring` and returns the force exerted by this spring.

```@example component
function spring_force(spring)
    -spring.k .* scalarize(spring.dir) .* (spring.x - spring.l) ./ spring.x
end
```

To create our system, we will first create the components: a mass and a spring. This is done as follows:

```@example component
m = 1.0
xy = [1.0, -1.0]
k = 1e4
l = 1.0
center = [0.0, 0.0]
g = [0.0, -9.81]
@named mass = Mass(m = m, xy = xy)
@named spring = Spring(k = k, l = l)
```

We can now create the equations describing this system, by connecting `spring` to `mass` and a fixed point.

```@example component
eqs = [connect_spring(spring, mass.pos, center)
       scalarize(D.(mass.v) .~ spring_force(spring) / mass.m .+ g)]
```

Finally, we can build the model using these equations and components.

```@example component
@named _model = ODESystem(eqs, t, [spring.x; spring.dir; mass.pos], [])
@named model = compose(_model, mass, spring)
```

We can take a look at the equations in the model using the `equations` function.

```@example component
equations(model)
```

The unknowns of this model are:

```@example component
unknowns(model)
```

And the parameters of this model are:

```@example component
parameters(model)
```

### Simplifying and solving this system

This system can be solved directly as a DAE using [one of the DAE solvers from DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/solvers/dae_solve/). However, we can symbolically simplify the system first beforehand. Running `structural_simplify` eliminates unnecessary variables from the model to give the leanest numerical representation of the system.

```@example component
sys = structural_simplify(model)
equations(sys)
```

We are left with only 4 equations involving 4 unknown variables (`mass.pos[1]`,
`mass.pos[2]`, `mass.v[1]`, `mass.v[2]`). We can solve the system by converting
it to an `ODEProblem`. Some observed variables are not expanded by default. To
view the complete equations, one can do

```@example component
full_equations(sys)
```

This is done as follows:

```@example component
prob = ODEProblem(sys, [], (0.0, 3.0))
sol = solve(prob, Rosenbrock23())
plot(sol)
```

What if we want the timeseries of a different variable? That information is not lost! Instead, `structural_simplify` simply changes unknown variables into `observed` variables.

```@example component
observed(sys)
```

These are explicit algebraic equations which can be used to reconstruct the required variables on the fly. This leads to dramatic computational savings since implicitly solving an ODE scales as O(n^3), so fewer unknowns are significantly better!

We can access these variables using the solution object. For example, let's retrieve the x-position of the mass over time:

```@example component
sol[mass.pos[1]]
```

We can also plot the path of the mass:

```@example component
plot(sol, idxs = (mass.pos[1], mass.pos[2]))
```
