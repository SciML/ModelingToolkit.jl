# Interactive simulation with a long-lived integrator

Interactive applications often receive a new input, advance the model for a short interval,
and update a plot. Initialize one integrator for the full simulation horizon and keep advancing
that integrator. Reinitializing between intervals is unnecessary and can rerun ModelingToolkit
initialization or discard solver state.

## Declaring the input

Model the external signal as an input **variable** and name it to `mtkcompile`, which
converts it into a parameter of the compiled system:

```julia
using ModelingToolkit
using OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D

@variables x(t) u(t)
sys = mtkcompile(System([D(x) ~ -x + u], t; name = :sys); inputs = [u])

is_parameter(sys, u)   # true -- `u` is now a parameter, driven from outside
```

This is the same `inputs` keyword used throughout the
[input-output](@ref inputoutput) and control tutorials, so a model written for
linearization, control design or `generate_control_function` needs no rewriting to be
stepped interactively. It matters most for models built from component libraries: there
the input arrives as an unconnected `RealInput` connector, and there is no top-level
`@parameters` declaration to reach for. Naming that connector variable to `mtkcompile`
turns it into a parameter you can assign between steps.

Without it, compiling a model whose input connector is left unconnected fails, because
the input is counted as an unknown with no equation to determine it:

```
ExtraVariablesSystemException: The system is unbalanced.
There are 8 highest order derivative variables and 7 equations.
```

If the signal is genuinely a scalar you declare yourself, `@parameters` works directly and
needs no `inputs` keyword:

```julia
using ModelingToolkit
using OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D

@variables x(t)
@parameters input
@mtkcompile sys = System([D(x) ~ -x + input], t)

prob = ODEProblem(sys, [x => 0.0, input => 0.0], (0.0, 4.0))
integrator = init(prob, Tsit5(); saveat = 0.0:0.05:4.0)

display_times = Float64[]
display_values = Float64[]
display_inputs = Float64[]
for level in (0.0, 1.0, -1.0, 0.0)
    integrator.ps[input] = level
    step!(integrator, 1.0, true)
    push!(display_times, integrator.t)
    push!(display_values, integrator[x])
    push!(display_inputs, level)
end

sol = integrator.sol
```

`integrator.ps[input] = level` uses the
[`SymbolicIndexingInterface`](https://docs.sciml.ai/SymbolicIndexingInterface/stable/)
parameter interface. It works with the parameter storage produced by structural
transformation and tells the integrator that its derivative data must be refreshed before
the next step. Avoid mutating the internal `integrator.p` storage directly.

`step!(integrator, interval, true)` advances by `interval` and stops exactly at its endpoint.
The integrator's solution accumulates the saved state history, while `integrator.t` and
`integrator[x]` provide the values needed to update live plot observables after each block.
Use a time span covering the intended session and express `saveat` in absolute simulation
times.

Use `reinit!` when the simulation really is restarting from a new initial condition or time,
not to divide one continuous run into display intervals.

## Streaming an input into a component-library model

The same keyword is what makes this work for models assembled from
[ModelingToolkitStandardLibrary](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/)
components. Leave the driving component's input connector unconnected and name it to
`mtkcompile`; the loop that feeds it is then identical to the scalar case.

```julia
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, unbound_inputs
using ModelingToolkitStandardLibrary.Mechanical.Translational: Mass, Spring, Damper, Fixed, Force

@named mass = Mass(m = 1)
@named spring = Spring(k = 100)
@named damper = Damper(d = 1)
@named fix = Fixed()
@named force = Force()          # `force.f` is left unconnected: it is the external input

eqs = [connect(fix.flange, spring.flange_a)
       connect(spring.flange_b, mass.flange)
       connect(fix.flange, damper.flange_a)
       connect(damper.flange_b, mass.flange)
       connect(mass.flange, force.flange)]
@named model = System(eqs, t; systems = [mass, spring, damper, fix, force])

sys = mtkcompile(model; inputs = [force.f.u])

prob = ODEProblem(sys, [mass.s => 0.0, mass.v => 0.0, force.f.u => 0.0], (0.0, 4.0))
integrator = init(prob, Tsit5(); saveat = 0.0:0.05:4.0)

for level in (0.0, 10.0, -10.0, 0.0)
    integrator.ps[force.f.u] = level
    step!(integrator, 1.0, true)
    # Update the display from the current symbolic values.
    current_position = integrator[mass.s]
    current_velocity = integrator[mass.v]
end
```

When you do not already know which variable carries the external signal — a model built by
someone else, or one assembled programmatically — [`unbound_inputs`](@ref ModelingToolkit.unbound_inputs) reports the input
variables that the connection structure leaves external, and its result can be passed
straight to `mtkcompile`:

```julia
unbound_inputs(model)                              # [force₊f₊u(t)]
sys = mtkcompile(model; inputs = unbound_inputs(model))
```

This inspects the connection graph of the uncompiled hierarchy, so call it before
`mtkcompile`. It is a heuristic rather than a guarantee: check what it returns, and name the
input explicitly when you already know it.

## Reading observed quantities from the integrator

A live display usually needs quantities that structural simplification eliminated from the
state vector. Those are available from the integrator through the same symbolic indexing
that a solution supports, so the display code does not have to reconstruct them:

```julia
@variables pos(t) vel(t) energy(t)
@parameters input
@mtkcompile osc = System(
    [D(pos) ~ vel, D(vel) ~ -pos + input, energy ~ (pos^2 + vel^2) / 2], t
)

oscprob = ODEProblem(osc, [pos => 1.0, vel => 0.0, input => 0.0], (0.0, 10.0))
integrator = init(oscprob, Tsit5())
step!(integrator, 1.0, true)

integrator[energy]                       # current value of an observed equation
integrator(integrator.t; idxs = energy)  # the same quantity from the step's interpolant
integrator([integrator.tprev, integrator.t]; idxs = [pos, energy])
```

`plot(integrator)` accepts the same `idxs` specifications as `plot(sol)`, including observed
equations, lists, phase-plane tuples and a plot function:

```julia
using Plots

plot(integrator; idxs = energy)
plot(integrator; idxs = [pos, energy])
plot(integrator; idxs = (pos, vel))
```

`plot(integrator)` draws only the step the integrator is currently on, which is what an
animation loop wants. Plot `integrator.sol` instead to draw the whole accumulated history.

## Parameters and time-varying inputs

A value declared with `@parameters` is treated as time-invariant when a saved solution evaluates
parameter-dependent observables. Mutating one during an interactive run updates future
dynamics, but it does not create a history of its previous values. Record the input in the
application, as the example does for display values, when only the live interface needs that
history.

When the input is part of the model's saved history, declare it with `@discretes input(t)`.
Discrete variables represent piecewise-constant, time-dependent values. Updates performed by
ModelingToolkit callbacks can be saved and later retrieved through symbolic solution indexing;
see [Saving discrete values](@ref save_discretes).

If all input-change times are known before the solve, a symbolic discrete event is usually a
better model than an external stepping loop. The long-lived-integrator pattern is intended for
inputs that arrive while the simulation is running, such as UI controls or measurements from
another process.
