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

## [(Experimental) Imperative affects](@id imp_affects)

The `ImperativeAffect` can be used as an alternative to the aforementioned functional affect form. Note
that `ImperativeAffect` is still experimental; to emphasize this, we do not export it and it should be
included as `ModelingToolkit.ImperativeAffect`. `ImperativeAffect` aims to simplify the manipulation of
system state.

We will use two examples to describe `ImperativeAffect`: a simple heater and a quadrature encoder.
These examples will also demonstrate advanced usage of `ModelingToolkit.SymbolicContinuousCallback`,
the low-level interface of the tuple form converts into that allows control over the SciMLBase-level
event that is generated for a continuous event.

### [Heater](@id heater_events)

Bang-bang control of a heater connected to a leaky plant requires hysteresis in order to prevent rapid control oscillation.

```@example events
@variables temp(t)
params = @parameters furnace_on_threshold=0.5 furnace_off_threshold=0.7 furnace_power=1.0 leakage=0.1 furnace_on(t)::Bool=false
eqs = [
    D(temp) ~ furnace_on * furnace_power - temp^2 * leakage
]
```

Our plant is simple. We have a heater that's turned on and off by the time-indexed parameter `furnace_on`
which adds `furnace_power` forcing to the system when enabled. We then leak heat proportional to `leakage`
as a function of the square of the current temperature.

We need a controller with hysteresis to control the plant. We wish the furnace to turn on when the temperature
is below `furnace_on_threshold` and off when above `furnace_off_threshold`, while maintaining its current state
in between. To do this, we create two continuous callbacks:

```@example events
using Setfield
furnace_disable = ModelingToolkit.SymbolicContinuousCallback(
    [temp ~ furnace_off_threshold],
    ModelingToolkit.ImperativeAffect(modified = (; furnace_on)) do x, o, c, i
        @set! x.furnace_on = false
    end)
furnace_enable = ModelingToolkit.SymbolicContinuousCallback(
    [temp ~ furnace_on_threshold],
    ModelingToolkit.ImperativeAffect(modified = (; furnace_on)) do x, o, c, i
        @set! x.furnace_on = true
    end)
```

We're using the explicit form of `SymbolicContinuousCallback` here, though
so far we aren't using anything that's not possible with the implicit interface.
You can also write

```julia
[temp ~ furnace_off_threshold] => ModelingToolkit.ImperativeAffect(modified = (;
    furnace_on)) do x, o, i, c
    @set! x.furnace_on = false
end
```

and it would work the same.

The `ImperativeAffect` is the larger change in this example. `ImperativeAffect` has the constructor signature

```julia
ImperativeAffect(f::Function; modified::NamedTuple, observed::NamedTuple, ctx)
```

that accepts the function to call, a named tuple of both the names of and symbolic values representing
values in the system to be modified, a named tuple of the values that are merely observed (that is, used from
the system but not modified), and a context that's passed to the affect function.

In our example, each event merely changes whether the furnace is on or off. Accordingly, we pass a `modified` tuple
`(; furnace_on)` (creating a `NamedTuple` equivalent to `(furnace_on = furnace_on)`). `ImperativeAffect` will then
evaluate this before calling our function to fill out all of the numerical values, then apply them back to the system
once our affect function returns. Furthermore, it will check that it is possible to do this assignment.

The function given to `ImperativeAffect` needs to have the signature:

```julia
f(modified::NamedTuple, observed::NamedTuple, ctx, integrator)::NamedTuple
```

The function `f` will be called with `observed` and `modified` `NamedTuple`s that are derived from their respective `NamedTuple` definitions.
In our example, if `furnace_on` is `false`, then the value of the `x` that's passed in as `modified` will be `(furnace_on = false)`.
The modified values should be passed out in the same format: to set `furnace_on` to `true` we need to return a tuple `(furnace_on = true)`.
The examples does this with Setfield, recreating the result tuple before returning it; the returned tuple may optionally be missing values as
well, in which case those values will not be written back to the problem.

Accordingly, we can now interpret the `ImperativeAffect` definitions to mean that when `temp = furnace_off_threshold` we
will write `furnace_on = false` back to the system, and when `temp = furnace_on_threshold` we will write `furnace_on = true` back
to the system.

```@example events
@named sys = ODESystem(
    eqs, t, [temp], params; continuous_events = [furnace_disable, furnace_enable])
ss = structural_simplify(sys)
prob = ODEProblem(ss, [temp => 0.0, furnace_on => true], (0.0, 10.0))
sol = solve(prob, Tsit5())
plot(sol)
hline!([sol.ps[furnace_off_threshold], sol.ps[furnace_on_threshold]],
    l = (:black, 1), primary = false)
```

Here we see exactly the desired hysteresis. The heater starts on until the temperature hits
`furnace_off_threshold`. The temperature then bleeds away until `furnace_on_threshold` at which
point the furnace turns on again until `furnace_off_threshold` and so on and so forth. The controller
is effectively regulating the temperature of the plant.

### [Quadrature Encoder](@id quadrature)

For a more complex application we'll look at modeling a quadrature encoder attached to a shaft spinning at a constant speed.
Traditionally, a quadrature encoder is built out of a code wheel that interrupts the sensors at constant intervals and two sensors slightly out of phase with one another.
A state machine can take the pattern of pulses produced by the two sensors and determine the number of steps that the shaft has spun. The state machine takes the new value
from each sensor and the old values and decodes them into the direction that the wheel has spun in this step.

```@example events
@variables theta(t) omega(t)
params = @parameters qA=0 qB=0 hA=0 hB=0 cnt::Int=0
eqs = [D(theta) ~ omega
       omega ~ 1.0]
```

Our continuous-time system is extremely simple. We have two unknown variables `theta` for the angle of the shaft
and `omega` for the rate at which it's spinning. We then have parameters for the state machine `qA, qB, hA, hB`
(corresponding to the current quadrature of the A/B sensors and the historical ones) and a step count `cnt`.

We'll then implement the decoder as a simple Julia function.

```@example events
function decoder(oldA, oldB, newA, newB)
    state = (oldA, oldB, newA, newB)
    if state == (0, 0, 1, 0) || state == (1, 0, 1, 1) || state == (1, 1, 0, 1) ||
       state == (0, 1, 0, 0)
        return 1
    elseif state == (0, 0, 0, 1) || state == (0, 1, 1, 1) || state == (1, 1, 1, 0) ||
           state == (1, 0, 0, 0)
        return -1
    elseif state == (0, 0, 0, 0) || state == (0, 1, 0, 1) || state == (1, 0, 1, 0) ||
           state == (1, 1, 1, 1)
        return 0
    else
        return 0 # err is interpreted as no movement
    end
end
```

Based on the current and old state, this function will return 1 if the wheel spun in the positive direction,
-1 if in the negative, and 0 otherwise.

The encoder state advances when the occlusion begins or ends. We model the
code wheel as simply detecting when `cos(100*theta)` is 0; if we're at a positive
edge of the 0 crossing, then we interpret that as occlusion (so the discrete `qA` goes to 1). Otherwise, if `cos` is
going negative, we interpret that as lack of occlusion (so the discrete goes to 0). The decoder function is
then invoked to update the count with this new information.

We can implement this in one of two ways: using edge sign detection or right root finding. For exposition, we
will implement each sensor differently.

For sensor A, we're using the edge detection method. By providing a different affect to `SymbolicContinuousCallback`'s
`affect_neg` argument, we can specify different behaviour for the negative crossing vs. the positive crossing of the root.
In our encoder, we interpret this as occlusion or nonocclusion of the sensor, update the internal state, and tick the decoder.

```@example events
qAevt = ModelingToolkit.SymbolicContinuousCallback([cos(100 * theta) ~ 0],
    ModelingToolkit.ImperativeAffect((; qA, hA, hB, cnt), (; qB)) do x, o, c, i
        @set! x.hA = x.qA
        @set! x.hB = o.qB
        @set! x.qA = 1
        @set! x.cnt += decoder(x.hA, x.hB, x.qA, o.qB)
        x
    end,
    affect_neg = ModelingToolkit.ImperativeAffect(
        (; qA, hA, hB, cnt), (; qB)) do x, o, c, i
        @set! x.hA = x.qA
        @set! x.hB = o.qB
        @set! x.qA = 0
        @set! x.cnt += decoder(x.hA, x.hB, x.qA, o.qB)
        x
    end)
```

The other way we can implement a sensor is by changing the root find.
Normally, we use left root finding; the affect will be invoked instantaneously _before_
the root is crossed. This makes it trickier to figure out what the new state is.
Instead, we can use right root finding:

```@example events
qBevt = ModelingToolkit.SymbolicContinuousCallback([cos(100 * theta - π / 2) ~ 0],
    ModelingToolkit.ImperativeAffect((; qB, hA, hB, cnt), (; qA, theta)) do x, o, c, i
        @set! x.hA = o.qA
        @set! x.hB = x.qB
        @set! x.qB = clamp(sign(cos(100 * o.theta - π / 2)), 0.0, 1.0)
        @set! x.cnt += decoder(x.hA, x.hB, o.qA, x.qB)
        x
    end; rootfind = SciMLBase.RightRootFind)
```

Here, sensor B is located `π / 2` behind sensor A in angular space, so we're adjusting our
trigger function accordingly. We here ask for right root finding on the callback, so we know
that the value of said function will have the "new" sign rather than the old one. Thus, we can
determine the new state of the sensor from the sign of the indicator function evaluated at the
affect activation point, with -1 mapped to 0.

We can now simulate the encoder.

```@example events
@named sys = ODESystem(
    eqs, t, [theta, omega], params; continuous_events = [qAevt, qBevt])
ss = structural_simplify(sys)
prob = ODEProblem(ss, [theta => 0.0], (0.0, pi))
sol = solve(prob, Tsit5(); dtmax = 0.01)
sol.ps[cnt]
```

`cos(100*theta)` will have 200 crossings in the half rotation we've gone through, so the encoder would notionally count 200 steps.
Our encoder counts 198 steps (it loses one step to initialization and one step due to the final state falling squarely on an edge).
