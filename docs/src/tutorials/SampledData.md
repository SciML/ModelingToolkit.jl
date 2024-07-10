```@setup clocks
using ModelingToolkit
using ModelingToolkit: t_nounits as t
```

# Clocks and Sampled-Data Systems

A sampled-data system contains both continuous-time and discrete-time components, such as a continuous-time plant model and a discrete-time control system. ModelingToolkit supports the modeling and simulation of sampled-data systems by means of *clocks*.

!!! danger "Experimental"
    
    The sampled-data interface is currently experimental and at any time subject to breaking changes **not** respecting semantic versioning.

!!! note "Negative shifts"
    
    The initial release of the sampled-data interface **only supports negative shifts**.

A clock can be seen as an *event source*, i.e., when the clock ticks, an event is generated. In response to the event the discrete-time logic is executed, for example, a control signal is computed. For basic modeling of sampled-data systems, the user does not have to interact with clocks explicitly, instead, the modeling is performed using the operators

  - [`Sample`](@ref)
  - [`Hold`](@ref)
  - [`ShiftIndex`](@ref)
  - [`SampleTime`](@ref)

When a continuous-time variable `x` is sampled using `xd = Sample(x, dt)`, the result is a discrete-time variable `xd` that is defined and updated whenever the clock ticks. `xd` is *only defined when the clock ticks*, which it does with an interval of `dt`. If `dt` is unspecified, the tick rate of the clock associated with `xd` is inferred from the context in which `xd` appears. Any variable taking part in the same equation as `xd` is inferred to belong to the same *discrete partition* as `xd`, i.e., belonging to the same clock. A system may contain multiple different discrete-time partitions, each with a unique clock. This allows for modeling of multi-rate systems and discrete-time processes located on different computers etc.

To make a discrete-time variable available to the continuous partition, the [`Hold`](@ref) operator is used. `xc = Hold(xd)` creates a continuous-time variable `xc` that is updated whenever the clock associated with `xd` ticks, and holds its value constant between ticks.

The operators [`Sample`](@ref) and [`Hold`](@ref) are thus providing the interface between continuous and discrete partitions.

The [`ShiftIndex`](@ref) operator is used to refer to past and future values of discrete-time variables. The example below illustrates its use, implementing the discrete-time system

```math
x(k+1) = 0.5x(k) + u(k)
y(k) = x(k)
```

```@example clocks
using ModelingToolkit
using ModelingToolkit: t_nounits as t
@variables x(t) y(t) u(t)
dt = 0.1                # Sample interval
clock = Clock(t, dt)    # A periodic clock with tick rate dt
k = ShiftIndex(clock)   # If the clock is not specified, clock inference will be used to try to infer the correct block from the equation context

eqs = [
    x(k) ~ 0.5x(k - 1) + u(k - 1),
    y ~ x
]
```

A few things to note in this basic example:

  - The `clock` is explicitly created and passed to `ShiftIndex`. When defining components, it may be beneficial to use a `ShiftIndex` without a clock to let the clock be inferred from context instead. The clock must be specified in at least one place in each clock partition.
  - The equation `x(k+1) = 0.5x(k) + u(k)` has been rewritten in terms of negative shifts since positive shifts are not yet supported.
  - `x` and `u` are automatically inferred to be discrete-time variables, since they appear in an equation with a discrete-time [`ShiftIndex`](@ref) `k`.
  - `y` is also automatically inferred to be a discrete-time-time variable, since it appears in an equation with another discrete-time variable `x`. `x,u,y` all belong to the same discrete-time partition, i.e., they are all updated at the same *instantaneous point in time* at which the clock ticks.
  - The equation `y ~ x` does not use any shift index, this is equivalent to `y(k) ~ x(k)`, i.e., discrete-time variables without shift index are assumed to refer to the variable at the current time step.
  - The equation `x(k) ~ 0.5x(k-1) + u(k-1)` indicates how `x` is updated, i.e., what the value of `x` will be at the *current* time step in terms of the *past* value. The output `y`, is given by the value of `x` at the *current* time step, i.e., `y(k) ~ x(k)`. If this logic was implemented in an imperative programming style, the logic would thus be

```julia
function discrete_step(x, u)
    x = 0.5x + u # x is updated to a new value, i.e., x(k) is computed
    y = x # y is assigned the current value of x, y(k) = x(k)
    return x, y # The state x now refers to x at the current time step, x(k), and y equals x, y(k) = x(k)
end
```

An alternative and *equivalent* way of writing the same system is

```@example clocks
eqs = [
    x(k + 1) ~ 0.5x(k) + u(k),
    y(k) ~ x(k)
]
```

but the use of positive time shifts is not yet supported. Instead, we *shifted all indices* by `-1` above, resulting in exactly the same difference equations. However, the next system is *not equivalent* to the previous one:

```@example clocks
eqs = [
    x(k) ~ 0.5x(k - 1) + u(k),
    y ~ x
]
```

In this last example, `u(k)` refers to the input at the new time point `k`., this system is equivalent to

```
eqs = [
    x(k+1) ~ 0.5x(k) + u(k+1),
    y(k) ~ x(k)
]
```

## Higher-order shifts

The expression `x(k-1)` refers to the value of `x` at the *previous* clock tick. Similarly, `x(k-2)` refers to the value of `x` at the clock tick before that. In general, `x(k-n)` refers to the value of `x` at the `n`th clock tick before the current one. As an example, the Z-domain transfer function

```math
H(z) = \dfrac{b_2 z^2 + b_1 z + b_0}{a_2 z^2 + a_1 z + a_0}
```

may thus be modeled as

```@example clocks
@parameters a0 a1 a2 b0 b1 b2
@variables y(t) [description = "Output"] u(t) [description = "Input"]
k = ShiftIndex(Clock(t, dt))
eqs = [
    a2 * y(k) + a1 * y(k - 1) + a0 * y(k - 2) ~ b2 * u(k) + b1 * u(k - 1) + b0 * u(k - 2)
]
nothing # hide
```

(see also [ModelingToolkitStandardLibrary](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/) for a discrete-time transfer-function component.)

## Initial conditions

The initial condition of discrete-time variables is defined using the [`ShiftIndex`](@ref) operator, for example

```julia
ODEProblem(model, [x(k) => 1.0], (0.0, 10.0))
```

If higher-order shifts are present, the corresponding initial conditions must be specified, e.g., the presence of the equation

```julia
x(k) ~ x(k - 1) + x(k - 2)
```

requires specification of the initial condition for both `x(k-1)` and `x(k-2)`.

## Multiple clocks

Multi-rate systems are easy to model using multiple different clocks. The following set of equations is valid, and defines *two different discrete-time partitions*, each with its own clock:

```julia
@variables r y yd1 yd2 ud1 ud2
@parameters kp
dt1, dt2 = 1, 3
yd1 ~ Sample(t, dt1)(y) # This automatically creates a `Clock(t, dt1)`
ud1 ~ kp * (Sample(t, dt1)(r) - yd1)
yd2 ~ Sample(t, dt2)(y)
ud2 ~ kp * (Sample(t, dt2)(r) - yd2)
nothing # hide
```

`yd1` and `ud1` belong to the same clock which ticks with an interval of `dt1`, while `yd2` and `ud2` belong to a different clock which ticks with an interval of `dt2`. The two clocks are *not synchronized*, i.e., they are not *guaranteed* to tick at the same point in time, even if one tick interval is a rational multiple of the other. Mechanisms for synchronization of clocks are not yet implemented.

In the example above, the clocks were implicitly created by the constructor to the `Sample` operator. Clocks can also be explicitly created using the constructor of the clock, e.g.,

```@example clocks
@variables y(t)
clock = Clock(t, 0.1)
Sample(clock)(y)
```

## Accessing discrete-time variables in the solution

Not yet available, work in progress.

## Referring to the sample interval in equations

It is common to want to refer to the sample interval when defining discrete-time equivalents of continuous-time components, such as integrators and derivatives etc. To facilitate this, the operator [`SampleTime`](@ref) is available, see the docstring for an example, as well as the definitions of discrete-time components in [ModelingToolkitStandardLibrary.Blocks](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/blocks/) for typical usage.

In the presence of multiple different sample rates, the value of the operator is given by the clock associated with the equation in which it appears. In the example below, the first appearance of `SampleTime()` returns `dt1` while the second returns `dt2`.

```@example clocks
yd1 ~ Sample(t, dt1)(y)
ud1 ~ kp * (Sample(t, dt1)(r) - yd1) * SampleTime()
yd2 ~ Sample(t, dt2)(y)
ud2 ~ kp * (Sample(t, dt2)(r) - yd2) * SampleTime()
```

## A complete example

Below, we model a simple continuous first-order system called `plant` that is controlled using a discrete-time controller `controller`. The reference signal is filtered using a discrete-time filter `filt` before being fed to the controller.

```@example clocks
using ModelingToolkit, Plots, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t
using ModelingToolkit: D_nounits as D
dt = 0.5 # Sample interval
clock = Clock(t, dt)
k = ShiftIndex(clock)

@mtkmodel Plant begin
    @variables begin
        x(t) = 1
        u(t) = 0
        y(t) = 0
    end
    @equations begin
        D(x) ~ -x + u
        y ~ x
    end
end

@mtkmodel Filt begin # Reference filter
    @variables begin
        x(t) = 0
        u(t) = 0
        y(t) = 0
    end
    begin
        a = 1 / exp(SampleTime()) # SampleTime() returns the sample interval of the clock
    end
    @equations begin
        x(k) ~ a * x(k - 1) + (1 - a) * u(k)
        y ~ x
    end
end

@mtkmodel Controller begin
    @parameters begin
        kp, [description = "Proportional gain"]
        ki, [description = "Integral gain"]
    end
    @variables begin
        y(t)
        r(t)
        ud(t)
        yd(t)
        e(t)
        I(t) = 0
    end
    @equations begin
        yd ~ Sample(y)
        e ~ (r - yd)
        I(k) ~ I(k - 1) + ki * e(k)
        ud ~ kp * e + I
    end
end

@mtkmodel ClosedLoop begin
    @components begin
        f = Filt()
        c = Controller(kp = 1, ki = 0.5)
        p = Plant()
    end
    @variables begin
        r(t)
    end
    @equations begin
        r ~ (t >= 5) - (t >= 15)  # reference signal
        f.u ~ r             # reference to filter input
        f.y ~ c.r           # filtered reference to controller reference
        Hold(c.ud) ~ p.u    # controller output to plant input (zero-order-hold)
        p.y ~ c.y           # plant output to controller feedback
    end
end

@mtkbuild model = ClosedLoop()
prob = ODEProblem(model, [], (0, 30))
sol = solve(prob, Tsit5(), kwargshandle = KeywordArgSilent) # Silence keywords required for now, to be removed in the future
plot(sol, title = "Continuous-time plant with discrete-time control system",
    titlefontsize = 10)
```

See [ModelingToolkitStandardLibrary.Blocks](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/blocks/) for discrete-time components with several additional features compared to the primitive example components above.

## Current limitations

As mentioned above, this funcitonality is still work in progress. The following known limitations are currently effective

  - Discrete-time variables cannot be accessed in the solution object.
  - Systems including clocked partitions cannot be linearized.
  - Discrete-time variables can only be indexed causally, i.e., `y(k) ~ y(k-1)` is valid while `y(k+1) ~ y(k)` is not.
