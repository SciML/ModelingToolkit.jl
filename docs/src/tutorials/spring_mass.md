# Component-Based Modeling a Spring-Mass System

In this tutorial we will build a simple component-based model of a spring-mass system. A spring-mass system consists of one or more masses connected by springs. [Hooke's law](https://en.wikipedia.org/wiki/Hooke%27s_law) gives the force exerted by a spring when it is extended or compressed by a given distance. This specifies a differential-equation system where the acceleration of the masses is specified using the forces acting on them.

## Copy-Paste Example

```julia
using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra

@variables t
D = Differential(t)

function Mass(; name, m = 1.0, xy = [0., 0.], u = [0., 0.])
    ps = @parameters m=m
    sts = @variables pos[1:2](t)=xy v[1:2](t)=u
    eqs = collect(D.(pos) .~ v)
    ODESystem(eqs, t, [pos..., v...], ps; name)
end

function Spring(; name, k = 1e4, l = 1.)
    ps = @parameters k=k l=l
    @variables x(t), dir[1:2](t)
    ODESystem(Equation[], t, [x, dir...], ps; name)
end

function connect_spring(spring, a, b)
    [
        spring.x ~ norm(collect(a .- b))
        collect(spring.dir .~ collect(a .- b))
    ]
end

spring_force(spring) = -spring.k .* collect(spring.dir) .* (spring.x - spring.l)  ./ spring.x

m = 1.0
xy = [1., -1.]
k = 1e4
l = 1.
center = [0., 0.]
g = [0., -9.81]
@named mass = Mass(m=m, xy=xy)
@named spring = Spring(k=k, l=l)

eqs = [
    connect_spring(spring, mass.pos, center)
    collect(D.(mass.v) .~ spring_force(spring) / mass.m .+ g)
]

@named _model = ODESystem(eqs, t)
@named model = compose(_model, mass, spring)
sys = structural_simplify(model)

prob = ODEProblem(sys, [], (0., 3.))
sol = solve(prob, Rosenbrock23())
plot(sol)
```

## Explanation
### Building the components
For each component we use a Julia function that returns an `ODESystem`. At the top, we define the fundamental properties of a `Mass`: it has a mass `m`, a position `pos` and a velocity `v`. We also define that the velocity is the rate of change of position with respect to time.

```julia
function Mass(; name, m = 1.0, xy = [0., 0.], u = [0., 0.])
    ps = @parameters m=m
    sts = @variables pos[1:2](t)=xy v[1:2](t)=u
    eqs = collect(D.(pos) .~ v)
    ODESystem(eqs, t, [pos..., v...], ps; name)
end
```

Note that this is an incompletely specified `ODESystem`. It cannot be simulated on its own since the equations for the velocity `v[1:2](t)` are unknown. Notice the addition of a `name` keyword. This allows us to generate different masses with different names. A `Mass` can now be constructed as:

```julia
Mass(name = :mass1)
```

Or using the `@named` helper macro

```julia
@named mass1 = Mass()
```

Next we build the spring component. It is characterised by the spring constant `k` and the length `l` of the spring when no force is applied to it. The state of a spring is defined by its current length and direction.

```julia
function Spring(; name, k = 1e4, l = 1.)
    ps = @parameters k=k l=l
    @variables x(t), dir[1:2](t)
    ODESystem(Equation[], t, [x, dir...], ps; name)
end
```

We now define functions that help construct the equations for a mass-spring system. First, the `connect_spring` function connects a `spring` between two positions `a` and `b`. Note that `a` and `b` can be the `pos` of a `Mass`, or just a fixed position such as `[0., 0.]`.

```julia
function connect_spring(spring, a, b)
    [
        spring.x ~ norm(collect(a .- b))
        collect(spring.dir .~ collect(a .- b))
    ]
end
```

Lastly, we define the `spring_force` function that takes a `spring` and returns the force exerted by this spring.

```julia
spring_force(spring) = -spring.k .* collect(spring.dir) .* (spring.x - spring.l)  ./ spring.x
```

To create our system, we will first create the components: a mass and a spring. This is done as follows:

```julia
m = 1.0
xy = [1., -1.]
k = 1e4
l = 1.
center = [0., 0.]
g = [0., -9.81]
@named mass = Mass(m=m, xy=xy)
@named spring = Spring(k=k, l=l)
```

We can now create the equations describing this system, by connecting `spring` to `mass` and a fixed point.

```julia
eqs = [
    connect_spring(spring, mass.pos, center)
    collect(D.(mass.v) .~ spring_force(spring) / mass.m .+ g)
]
```

Finally, we can build the model using these equations and components.

```julia
@named _model = ODESystem(eqs, t)
@named model = compose(_model, mass, spring)
```

We can take a look at the equations in the model using the `equations` function.

```julia
equations(model)

7-element Vector{Equation}:
 Differential(t)(mass₊v[1](t)) ~ -spring₊k*spring₊dir[1](t)*(mass₊m^-1)*(spring₊x(t) - spring₊l)*(spring₊x(t)^-1)
 Differential(t)(mass₊v[2](t)) ~ -9.81 - (spring₊k*spring₊dir[2](t)*(mass₊m^-1)*(spring₊x(t) - spring₊l)*(spring₊x(t)^-1))
 spring₊x(t) ~ sqrt(abs2(mass₊pos[1](t)) + abs2(mass₊pos[2](t)))
 spring₊dir[1](t) ~ mass₊pos[1](t)
 spring₊dir[2](t) ~ mass₊pos[2](t)
 Differential(t)(mass₊pos[1](t)) ~ mass₊v[1](t)
 Differential(t)(mass₊pos[2](t)) ~ mass₊v[2](t)
```

The states of this model are:

```julia
states(model)

7-element Vector{Term{Real, Base.ImmutableDict{DataType, Any}}}:
 mass₊v[1](t)
 mass₊v[2](t)
 spring₊x(t)
 mass₊pos[1](t)
 mass₊pos[2](t)
 spring₊dir[1](t)
 spring₊dir[2](t)
```

And the parameters of this model are:

```julia
parameters(model)

6-element Vector{Sym{Real, Base.ImmutableDict{DataType, Any}}}:
 spring₊k
 mass₊m
 spring₊l
 mass₊m
 spring₊k
 spring₊l
```

### Simplifying and solving this system

This system can be solved directly as a DAE using [one of the DAE solvers from DifferentialEquations.jl](https://diffeq.sciml.ai/stable/solvers/dae_solve/). However, we can symbolically simplify the system first beforehand. Running `structural_simplify` eliminates unnecessary variables from the model to give the leanest numerical representation of the system.

```julia
sys = structural_simplify(model)
equations(sys)

4-element Vector{Equation}:
 Differential(t)(mass₊v[1](t)) ~ -spring₊k*mass₊pos[1](t)*(mass₊m^-1)*(sqrt(abs2(mass₊pos[1](t)) + abs2(mass₊pos[2](t))) - spring₊l)*(sqrt(abs2(mass₊pos[1](t)) + abs2(mass₊pos[2](t)))^-1)
 Differential(t)(mass₊v[2](t)) ~ -9.81 - (spring₊k*mass₊pos[2](t)*(mass₊m^-1)*(sqrt(abs2(mass₊pos[1](t)) + abs2(mass₊pos[2](t))) - spring₊l)*(sqrt(abs2(mass₊pos[1](t)) + abs2(mass₊pos[2](t)))^-1))
 Differential(t)(mass₊pos[1](t)) ~ mass₊v[1](t)
 Differential(t)(mass₊pos[2](t)) ~ mass₊v[2](t)
```

We are left with only 4 equations involving 4 state variables (`mass.pos[1]`, `mass.pos[2]`, `mass.v[1]`, `mass.v[2]`). We can solve the system by converting it to an `ODEProblem` in mass matrix form and solving with an [`ODEProblem` mass matrix solver](https://diffeq.sciml.ai/stable/solvers/dae_solve/#OrdinaryDiffEq.jl-(Mass-Matrix)). This is done as follows:

```julia
prob = ODEProblem(sys, [], (0., 3.))
sol = solve(prob, Rosenbrock23())
plot(sol)
```

What if we want the timeseries of a different variable? That information is not lost! Instead, `structural_simplify` simply changes state variables into `observed` variables.

```julia
observed(sys)

3-element Vector{Equation}:
 spring₊dir[2](t) ~ mass₊pos[2](t)
 spring₊dir[1](t) ~ mass₊pos[1](t)
 spring₊x(t) ~ sqrt(abs2(mass₊pos[1](t)) + abs2(mass₊pos[2](t)))
```

These are explicit algebraic equations which can be used to reconstruct the required variables on the fly. This leads to dramatic computational savings since implicitly solving an ODE scales as O(n^3), so fewer states are signficantly better!

We can access these variables using the solution object. For example, let's retrieve the length of the spring over time:

```julia
sol[spring.x]
```

We can also plot its timeseries:

```julia
plot(sol, vars = [spring.x])
```