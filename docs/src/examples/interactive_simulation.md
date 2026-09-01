# Interactive simulation with a long-lived integrator

Interactive applications often receive a new input, advance the model for a short interval,
and update a plot. Initialize one integrator for the full simulation horizon and keep advancing
that integrator. Reinitializing between intervals is unnecessary and can rerun ModelingToolkit
initialization or discard solver state.

The following model uses a parameter as an externally controlled input:

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

## Parameters and time-varying inputs

An ordinary `@parameter` is treated as time-invariant when a saved solution evaluates
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
