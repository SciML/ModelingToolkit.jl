# [Event Handling and Callback Functions](@id events)

ModelingToolkit provides several ways to represent system events, which enable
system state or parameters to be changed when certain conditions are satisfied,
or can be used to detect discontinuities. These events are ultimately converted
into DifferentialEquations.jl [`ContinuousCallback`s or
`DiscreteCallback`s](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/),
or into more specialized callback types from the
[DiffEqCallbacks.jl](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_library/)
library.

[`ODESystem`](@ref)s and [`SDESystem`](@ref)s accept keyword arguments
`continuous_events` and `discrete_events` to symbolically encode continuous or
discrete callbacks. [`JumpSystem`](@ref)s currently support only
`discrete_events`. Continuous events are applied when a given condition becomes
zero, with root finding used to determine the time at which a zero crossing
occurred. Discrete events are applied when a condition tested after each
timestep evaluates to true. See the [DifferentialEquations
docs](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/)
for more detail.

Events involve both a *condition* function (for the zero crossing or truth
test), and an *affect* function (for determining how to update the system when
the event occurs). These can both be specified symbolically, but a more [general
functional affect](@ref func_affects) representation is also allowed, as described
below.

## Continuous Events

The basic purely symbolic continuous event interface to encode *one* continuous
event is

```julia
AbstractSystem(eqs, ...; continuous_events::Vector{Equation})
AbstractSystem(eqs, ...; continuous_events::Pair{Vector{Equation}, Vector{Equation}})
```

In the former, equations that evaluate to 0 will represent conditions that should
be detected by the integrator, for example to force stepping to times of
discontinuities. The latter allow modeling of events that have an effect on the
state, where the first entry in the `Pair` is a vector of equations describing
event conditions, and the second vector of equations describes the effect on the
state. Each affect equation must be of the form

```julia
single_unknown_variable ~ expression_involving_any_variables_or_parameters
```

or

```julia
single_parameter ~ expression_involving_any_variables_or_parameters
```

In this basic interface, multiple variables can be changed in one event, or
multiple parameters, but not a mix of parameters and variables. The latter can
be handled via more [general functional affects](@ref func_affects).

Finally, multiple events can be encoded via a `Vector{Pair{Vector{Equation}, Vector{Equation}}}`.

### Example: Friction

The system below illustrates how continuous events can be used to model Coulomb
friction

```@example events
using ModelingToolkit, OrdinaryDiffEq, Plots
using ModelingToolkit: t_nounits as t, D_nounits as D

function UnitMassWithFriction(k; name)
    @variables x(t)=0 v(t)=0
    eqs = [D(x) ~ v
           D(v) ~ sin(t) - k * sign(v)]
    ODESystem(eqs, t; continuous_events = [v ~ 0], name) # when v = 0 there is a discontinuity
end
@mtkbuild m = UnitMassWithFriction(0.7)
prob = ODEProblem(m, Pair[], (0, 10pi))
sol = solve(prob, Tsit5())
plot(sol)
```

### Example: Bouncing ball

In the documentation for
[DifferentialEquations](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/#Example-1:-Bouncing-Ball),
we have an example where a bouncing ball is simulated using callbacks which have
an `affect!` on the state. We can model the same system using ModelingToolkit
like this

```@example events
@variables x(t)=1 v(t)=0

root_eqs = [x ~ 0]  # the event happens at the ground x(t) = 0
affect = [v ~ -v] # the effect is that the velocity changes sign

@mtkbuild ball = ODESystem([D(x) ~ v
                            D(v) ~ -9.8], t; continuous_events = root_eqs => affect) # equation => affect

tspan = (0.0, 5.0)
prob = ODEProblem(ball, Pair[], tspan)
sol = solve(prob, Tsit5())
@assert 0 <= minimum(sol[x]) <= 1e-10 # the ball never went through the floor but got very close
plot(sol)
```

### Test bouncing ball in 2D with walls

Multiple events? No problem! This example models a bouncing ball in 2D that is enclosed by two walls at $y = \pm 1.5$.

```@example events
@variables x(t)=1 y(t)=0 vx(t)=0 vy(t)=2

continuous_events = [[x ~ 0] => [vx ~ -vx]
                     [y ~ -1.5, y ~ 1.5] => [vy ~ -vy]]

@mtkbuild ball = ODESystem(
    [
        D(x) ~ vx,
        D(y) ~ vy,
        D(vx) ~ -9.8 - 0.1vx, # gravity + some small air resistance
        D(vy) ~ -0.1vy
    ], t; continuous_events)

tspan = (0.0, 10.0)
prob = ODEProblem(ball, Pair[], tspan)

sol = solve(prob, Tsit5())
@assert 0 <= minimum(sol[x]) <= 1e-10 # the ball never went through the floor but got very close
@assert minimum(sol[y]) > -1.5 # check wall conditions
@assert maximum(sol[y]) < 1.5  # check wall conditions

tv = sort([LinRange(0, 10, 200); sol.t])
plot(sol(tv)[y], sol(tv)[x], line_z = tv)
vline!([-1.5, 1.5], l = (:black, 5), primary = false)
hline!([0], l = (:black, 5), primary = false)
```

### [Generalized functional affect support](@id func_affects)

In some instances, a more flexible response to events is needed, which cannot be
encapsulated by symbolic equations. For example, a component may implement
complex behavior that is inconvenient or impossible to represent symbolically.
ModelingToolkit therefore supports regular Julia functions as affects: instead
of one or more equations, an affect is defined as a `tuple`:

```julia
[x ~ 0] => (affect!, [v, x], [p, q], [discretes...], ctx)
```

where, `affect!` is a Julia function with the signature: `affect!(integ, u, p, ctx)`; `[u,v]` and `[p,q]` are the symbolic unknowns (variables) and parameters
that are accessed by `affect!`, respectively; `discretes` are the parameters modified by `affect!`, if any;
and `ctx` is any context that is passed to `affect!` as the `ctx` argument.

`affect!` receives a [DifferentialEquations.jl
integrator](https://docs.sciml.ai/DiffEqDocs/stable/basics/integrator/)
as its first argument, which can then be used to access unknowns and parameters
that are provided in the `u` and `p` arguments (implemented as `NamedTuple`s).
The integrator can also be manipulated more generally to control solution
behavior, see the [integrator
interface](https://docs.sciml.ai/DiffEqDocs/stable/basics/integrator/)
documentation. In affect functions, we have that

```julia
function affect!(integ, u, p, ctx)
    # integ.t is the current time
    # integ.u[u.v] is the value of the unknown `v` above
    # integ.ps[p.q] is the value of the parameter `q` above
end
```

When accessing variables of a sub-system, it can be useful to rename them
(alternatively, an affect function may be reused in different contexts):

```julia
[x ~ 0] => (affect!, [resistor₊v => :v, x], [p, q => :p2], [], ctx)
```

Here, the symbolic variable `resistor₊v` is passed as `v` while the symbolic
parameter `q` has been renamed `p2`.

As an example, here is the bouncing ball example from above using the functional
affect interface:

```@example events
sts = @variables x(t), v(t)
par = @parameters g = 9.8
bb_eqs = [D(x) ~ v
          D(v) ~ -g]

function bb_affect!(integ, u, p, ctx)
    integ.u[u.v] = -integ.u[u.v]
end

reflect = [x ~ 0] => (bb_affect!, [v], [], [], nothing)

@mtkbuild bb_sys = ODESystem(bb_eqs, t, sts, par,
    continuous_events = reflect)

u0 = [v => 0.0, x => 1.0]

bb_prob = ODEProblem(bb_sys, u0, (0, 5.0))
bb_sol = solve(bb_prob, Tsit5())

plot(bb_sol)
```

## Discrete events support

In addition to continuous events, discrete events are also supported. The
general interface to represent a collection of discrete events is

```julia
AbstractSystem(eqs, ...; discrete_events = [condition1 => affect1, condition2 => affect2])
```

where conditions are symbolic expressions that should evaluate to `true` when an
individual affect should be executed. Here `affect1` and `affect2` are each
either a vector of one or more symbolic equations, or a functional affect, just
as for continuous events. As before, for any *one* event the symbolic affect
equations can either all change unknowns (i.e. variables) or all change
parameters, but one cannot currently mix unknown variable and parameter changes within one
individual event.

### Example: Injecting cells into a population

Suppose we have a population of `N(t)` cells that can grow and die, and at time
`t1` we want to inject `M` more cells into the population. We can model this by

```@example events
@parameters M tinject α
@variables N(t)
Dₜ = Differential(t)
eqs = [Dₜ(N) ~ α - N]

# at time tinject we inject M cells
injection = (t == tinject) => [N ~ N + M]

u0 = [N => 0.0]
tspan = (0.0, 20.0)
p = [α => 100.0, tinject => 10.0, M => 50]
@mtkbuild osys = ODESystem(eqs, t, [N], [α, M, tinject]; discrete_events = injection)
oprob = ODEProblem(osys, u0, tspan, p)
sol = solve(oprob, Tsit5(); tstops = 10.0)
plot(sol)
```

Notice, with generic discrete events that we want to occur at one or more fixed
times, we need to also set the `tstops` keyword argument to `solve` to ensure
the integrator stops at that time. In the next section, we show how one can
avoid this by using a preset-time callback.

Note that more general logical expressions can be built, for example, suppose we
want the event to occur at that time only if the solution is smaller than 50% of
its steady-state value (which is 100). We can encode this by modifying the event
to

```@example events
injection = ((t == tinject) & (N < 50)) => [N ~ N + M]

@mtkbuild osys = ODESystem(eqs, t, [N], [M, tinject, α]; discrete_events = injection)
oprob = ODEProblem(osys, u0, tspan, p)
sol = solve(oprob, Tsit5(); tstops = 10.0)
plot(sol)
```

Since the solution is *not* smaller than half its steady-state value at the
event time, the event condition now returns false. Here we used logical and,
`&`, instead of the short-circuiting logical and, `&&`, as currently the latter
cannot be used within symbolic expressions.

Let's now also add a drug at time `tkill` that turns off production of new
cells, modeled by setting `α = 0.0`

```@example events
@parameters tkill

# we reset the first event to just occur at tinject
injection = (t == tinject) => [N ~ N + M]

# at time tkill we turn off production of cells
killing = (t == tkill) => [α ~ 0.0]

tspan = (0.0, 30.0)
p = [α => 100.0, tinject => 10.0, M => 50, tkill => 20.0]
@mtkbuild osys = ODESystem(eqs, t, [N], [α, M, tinject, tkill];
    discrete_events = [injection, killing])
oprob = ODEProblem(osys, u0, tspan, p)
sol = solve(oprob, Tsit5(); tstops = [10.0, 20.0])
plot(sol)
```

### Periodic and preset-time events

Two important subclasses of discrete events are periodic and preset-time
events.

A preset-time event is triggered at specific set times, which can be
passed in a vector like

```julia
discrete_events = [[1.0, 4.0] => [v ~ -v]]
```

This will change the sign of `v` *only* at `t = 1.0` and `t = 4.0`.

As such, our last example with treatment and killing could instead be modeled by

```@example events
injection = [10.0] => [N ~ N + M]
killing = [20.0] => [α ~ 0.0]

p = [α => 100.0, M => 50]
@mtkbuild osys = ODESystem(eqs, t, [N], [α, M];
    discrete_events = [injection, killing])
oprob = ODEProblem(osys, u0, tspan, p)
sol = solve(oprob, Tsit5())
plot(sol)
```

Notice, one advantage of using a preset-time event is that one does not need to
also specify `tstops` in the call to solve.

A periodic event is triggered at fixed intervals (e.g. every Δt seconds). To
specify a periodic interval, pass the interval as the condition for the event.
For example,

```julia
discrete_events = [1.0 => [v ~ -v]]
```

will change the sign of `v` at `t = 1.0`, `2.0`, ...

Finally, we note that to specify an event at precisely one time, say 2.0 below,
one must still use a vector

```julia
discrete_events = [[2.0] => [v ~ -v]]
```

## Saving discrete values

Time-dependent parameters which are updated in callbacks are termed as discrete variables.
ModelingToolkit enables automatically saving the timeseries of these discrete variables,
and indexing the solution object to obtain the saved timeseries. Consider the following
example:

```@example events
@variables x(t)
@parameters c(t)

@mtkbuild sys = ODESystem(
    D(x) ~ c * cos(x), t, [x], [c]; discrete_events = [1.0 => [c ~ c + 1]])

prob = ODEProblem(sys, [x => 0.0], (0.0, 2pi), [c => 1.0])
sol = solve(prob, Tsit5())
sol[c]
```

The solution object can also be interpolated with the discrete variables

```@example events
sol([1.0, 2.0], idxs = [c, c * cos(x)])
```

Note that only time-dependent parameters will be saved. If we repeat the above example with
this change:

```@example events
@variables x(t)
@parameters c

@mtkbuild sys = ODESystem(
    D(x) ~ c * cos(x), t, [x], [c]; discrete_events = [1.0 => [c ~ c + 1]])

prob = ODEProblem(sys, [x => 0.0], (0.0, 2pi), [c => 1.0])
sol = solve(prob, Tsit5())
sol.ps[c] # sol[c] will error, since `c` is not a timeseries value
```

It can be seen that the timeseries for `c` is not saved.
