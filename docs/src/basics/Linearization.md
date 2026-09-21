# [Linearization](@id linearization)

A nonlinear dynamical system with state (differential and algebraic) ``x`` and input signals ``u``

```math
M \dot x = f(x, u)
```

can be linearized using the function [`linearize`](@ref) to produce a linear statespace system on the form

```math
\begin{aligned}
\dot x &= Ax + Bu\\
y &= Cx + Du
\end{aligned}
```

The `linearize` function expects the user to specify the inputs ``u`` and the outputs ``y`` using the syntax shown in the example below. The system model is *not* supposed to be simplified before calling `linearize`:

## Example

```@example LINEARIZE
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
@variables x(t)=0 y(t) u(t) r(t)=0
@parameters kp = 1

eqs = [u ~ kp * (r - y) # P controller
       D(x) ~ -x + u    # First-order plant
       y ~ x]           # Output equation

@named sys = System(eqs, t) # Do not call @mtkcompile when linearizing
matrices, simplified_sys = linearize(sys, [r], [y]) # Linearize from r to y
matrices
```

The named tuple `matrices` contains the matrices of the linear statespace representation, while `simplified_sys` is an `System` that, among other things, indicates the unknown variable order in the linear system through

```@example LINEARIZE
using ModelingToolkit: inputs, outputs
[unknowns(simplified_sys); inputs(simplified_sys); outputs(simplified_sys)]
```

!!! note "Inputs must be unconnected"

    The model above has 4 variables but only three equations, there is no equation specifying the value of `r` since `r` is an input. This means that only unbalanced models can be linearized, or in other words, models that are balanced and can be simulated _cannot_ be linearized. To learn more about this, see [How to linearize a ModelingToolkit model (YouTube)](https://www.youtube.com/watch?v=-XOux-2XDGI&t=395s). Also see [ModelingToolkitStandardLibrary: Linear analysis](https://docs.sciml.ai/ModelingToolkit/stable/tutorials/linear_analysis/) for utilities that make linearization of completed models easier.

!!! note "Un-simplified system"

    Linearization expects `sys` to be un-simplified, i.e., `mtkcompile` or `@mtkcompile` should not be called on the system before linearizing.

## Operating point

The operating point to linearize around can be specified with the keyword argument `op` like this: `op = Dict(x => 1, r => 2)`. The operating point may include specification of unknown variables, input variables and parameters. For variables that are not specified in `op`, the default value specified in the model will be used if available, if no value is specified, an error is thrown.

## Batch linearization and algebraic variables

If linearization is to be performed around multiple operating points, the simplification of the system has to be carried out a single time only. To facilitate this, the lower-level function [`ModelingToolkit.linearization_function`](@ref) is available. This function further allows you to obtain separate Jacobians for the differential and algebraic parts of the model. For ODE models without algebraic equations, the statespace representation above is available from the output of `linearization_function` as `A, B, C, D = f_x, f_u, h_x, h_u`.

All variables that will be fixed by an operating point _must_ be provided in the operating point to `linearization_function`. For example, if the operating points fix the value of
`x`, `y` and `z` then an operating point with constant values for these variables (e.g. `Dict(x => 1.0, y => 1.0, z => 1.0)`) must be provided. The constant values themselves
do not matter and can be changed by subsequent operating points.

One approach to batch linearization would be to call `linearize` in a loop, providing a different operating point each time. For example:

```@example LINEARIZE
using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Blocks

@parameters k=10 k3=2 c=1
@variables x(t)=0 [bounds = (-0.5, 1.5)]
@variables v(t) = 0

@named y = Blocks.RealOutput()
@named u = Blocks.RealInput()

eqs = [D(x) ~ v
       D(v) ~ -k * x - k3 * x^3 - c * v + 10u.u
       y.u ~ x]

@named duffing = System(eqs, t, systems = [y, u], initial_conditions = [u.u => 0])

# pass a constant value for `x`, since it is the variable we will change in operating points
linfun, simplified_sys = linearization_function(duffing, [u.u], [y.u]; op = Dict(x => NaN));

println(linearize(simplified_sys, linfun; op = Dict(x => 1.0)))
println(linearize(simplified_sys, linfun; op = Dict(x => 0.0)))

@time linearize(simplified_sys, linfun; op = Dict(x => 0.0))

nothing # hide
```

However, this route is still expensive since it has to repeatedly process the symbolic map provided to `op`. `linearize` is simply a wrapper for creating and solving a
[`ModelingToolkit.LinearizationProblem`](@ref). This object is symbolically indexable, and can thus integrate with SymbolicIndexingInterface.jl for fast updates.

```@example LINEARIZE
using SymbolicIndexingInterface

# The second argument is the value of the independent variable `t`.
linprob = LinearizationProblem(linfun, 1.0)
# It can be mutated
linprob.t = 0.0
# create a setter function to update `x` efficiently
setter! = setu(linprob, x)

function fast_linearize!(problem, setter!, value)
    setter!(problem, value)
    solve(problem)
end

println(fast_linearize!(linprob, setter!, 1.0))
println(fast_linearize!(linprob, setter!, 0.0))

@time fast_linearize!(linprob, setter!, 1.0)

nothing # hide
```

Note that `linprob` above can be interacted with similar to a normal `ODEProblem`.

```@repl LINEARIZE
prob[x]
prob[x] = 1.5
prob[x]
```

## Symbolic linearization

The function [`ModelingToolkit.linearize_symbolic`](@ref) works similar to [`ModelingToolkit.linearize`](@ref) but returns symbolic rather than numeric Jacobians. Symbolic linearization have several limitations and no all systems that can be linearized numerically can be linearized symbolically.

## Input derivatives

Physical systems are always *proper*, i.e., they do not differentiate causal inputs. However, ModelingToolkit allows you to model non-proper systems, such as inverse models, and may sometimes fail to find a realization of a proper system on proper form. In these situations, `linearize` may throw an error mentioning

```
Input derivatives appeared in expressions (-g_z\g_u != 0)
```

This means that to simulate this system, some order of derivatives of the input is required. To allow `linearize` to proceed in this situation, one may pass the keyword argument `allow_input_derivatives = true`, in which case the resulting model will have twice as many inputs, ``2n_u``, where the last ``n_u`` inputs correspond to ``\dot u``.

If the modeled system is actually proper (but MTK failed to find a proper realization), further numerical simplification can be applied to the resulting statespace system to obtain a proper form. Such simplification is currently available in the package [ControlSystemsMTK](https://juliacontrol.github.io/ControlSystemsMTK.jl/dev/#Internals:-Transformation-of-non-proper-models-to-proper-statespace-form).

## Tools for linear analysis

ModelingToolkit contains a set of [tools for more advanced linear analysis](https://docs.sciml.ai/ModelingToolkit/stable/tutorials/linear_analysis/). These can be used to make it easier to work with and analyze causal models, such as control and signal-processing systems.

Also see [ControlSystemsMTK.jl](https://juliacontrol.github.io/ControlSystemsMTK.jl/dev/) for an interface to [ControlSystems.jl](https://github.com/JuliaControl/ControlSystems.jl) that contains tools for linear analysis and frequency-domain analysis.

## Hybrid systems with clock partitions

`linearize` accounts for the continuous equations of a model only. A model that contains
synchronous (clocked) subsystems, written with `ShiftIndex`, `Sample` and `Hold` or built
from components such as those in DiscreteComponents, is linearized with
[`linearize_hybrid`](@ref). The model is split into its clock partitions, one per periodic
clock plus the continuous partition, and every partition is linearized separately. Each signal
that crosses a clock boundary becomes an input of the partition it enters and an output of
the partition it leaves, in addition to the inputs and outputs requested by the user.

```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

dt = 0.1
@variables x(t) y(t) u(t) d(t) yd(t) ud(t)
@parameters kp = 2.0
eqs = [
    yd ~ Sample(dt)(y)      # sample the plant output
    ud ~ -kp * yd           # discrete proportional controller
    u ~ Hold(ud)            # hold the controller output
    D(x) ~ -x + u + d       # plant with disturbance input d
    y ~ x
]
@named sys = System(eqs, t)
hl = linearize_hybrid(sys, [d], [y]; op = Dict(x => 0.0, d => 0.0))
```

The result is a [`HybridLinearization`](@ref). The continuous partition comes first, followed
by the discrete partitions.

```julia
cont, disc = hl.partitions
cont.inputs   # [d, Hold(ud)]: the disturbance and the held controller output
cont.outputs  # [y]: the plant output, which is also sampled by the controller
disc.inputs   # [Sample(dt)(y)]
disc.outputs  # [ud]
disc.Ts       # 0.1
hl.connections
```

The `inputs` and `outputs` of a partition list the user-specified signals first, in the order
given by the user; the fields `nu_user` and `ny_user` give the number of user-specified
inputs and outputs. The remaining signals are the boundary signals, grouped by the partition
they come from or go to. The connections state which output of which partition drives which
input of which other partition, `(; from, output, to, input)`. [`clock_boundary`](@ref)
collects the connections from one partition to another as index vectors, `outputs` into the
outputs of the partition the signals leave and `inputs` into the inputs of the partition they
enter, which are the index sets required by the advanced interface of
`ControlSystemsBase.feedback`. The variables of the boundary signals can be used as signal
names with `RobustAndOptimalControl.connect`.

For a discrete partition with sample interval `Ts`, the state consists of the values of the
discrete variables at the previous tick, denoted `xₜ₋₁` in `unknowns`, and the matrices
describe the update performed at a tick,

```math
\begin{aligned}
x(k) &= A x(k-1) + B u(k) \\
y(k) &= C x(k-1) + D u(k)
\end{aligned}
```

whose transfer function is ``C(zI - A)^{-1}B + D``. The continuous partition has the usual
form ``\dot x = Ax + Bu``, ``y = Cx + Du``.

To analyze the sampled-data loop, the partitions are brought to a common time domain and
interconnected according to `hl.connections`. The example below discretizes the continuous
partition with zero-order hold, which is the exact description of the loop at the sampling
instants when all discrete partitions share one clock, and closes the loop with
`ControlSystemsBase.feedback`.

```julia
using ControlSystemsBase
Pd = c2d(ss(cont.A, cont.B, cont.C, cont.D), disc.Ts)
Cd = ss(disc.A, disc.B, disc.C, disc.D, disc.Ts)
plant_to_ctrl = clock_boundary(hl, 1, 2)  # the plant output, sampled by the controller
ctrl_to_plant = clock_boundary(hl, 2, 1)  # the controller output, held and applied to the plant
# The boundary signals close the loop; the user-specified signals of the plant remain external.
G = feedback(
    Pd, Cd; Y1 = plant_to_ctrl.outputs, U2 = plant_to_ctrl.inputs,
    Y2 = ctrl_to_plant.outputs, U1 = ctrl_to_plant.inputs,
    W1 = 1:cont.nu_user, Z1 = 1:cont.ny_user, pos_feedback = true
)
```

This interconnection is the Redheffer star product of the two partitions over the boundary
signals. The functions `starprod` and `lft` of ControlSystemsBase compute it when the
interconnected signals are the last inputs and outputs of the first system and the first
inputs and outputs of the second, and when every output is either external or
interconnected. The partitions list the user-specified signals first, and a signal may be
both a user-specified output and a boundary output, as `y` is here, so the index sets are
passed to `feedback` explicitly in general. In the example above, `lft(Pd, Cd)` is not
applicable because `y` is the only output of the plant, whereas a plant whose sampled output
differs from its user-specified output is assembled with `lft(Pd, Cd)` directly when the
controller has no user-specified signals.

The discrete-time model describes the loop at the sampling instants. The content of the
sampled signals above the Nyquist frequency of the sampler, ``\omega_N = \pi / T_s``, is
folded by the sampling and is not represented by the model. Whether such content is present
can be checked on the partition that produces the sampled signals. Every signal entering a
sampler is an output of the continuous partition, driven by its inputs: the external inputs
and the held signals, whose piecewise constant form carries images of the discrete signal at
multiples of the sampling frequency. If the gain from all inputs of the continuous partition
to the sampled outputs is small above ``\omega_N`` compared with the gain below it, the
sampled signals carry little content that the sampling folds, and the discrete-time model of
the loop is representative.

```julia
P = ss(cont.A, cont.B, cont.C, cont.D)
Gb = P[plant_to_ctrl.outputs, :]  # from all inputs of the plant to the sampled signals
ωN = pi / disc.Ts
ω_above = exp10.(range(log10(ωN), log10(1000ωN), length = 200))
ω_below = exp10.(range(log10(ωN) - 3, log10(ωN), length = 200))
gain_above = maximum(abs, freqresp(Gb, ω_above))
gain_below = maximum(abs, freqresp(Gb, ω_below))
gain_above / gain_below
```

For this plant the ratio is about 0.03, that is, the plant attenuates the content above the
Nyquist frequency by a factor of thirty relative to its passband, and the sampled plant output
is well described by the discrete-time model. For partitions with several sampled signals or
inputs, the singular values (`sigma`) are the corresponding measure.

In a multi-rate model, a signal leaving a discrete partition for a partition on a slower
clock is decimated at the slower clock, which folds the band between the two Nyquist
frequencies, ``\pi / T_{s,\mathrm{slow}}`` to ``\pi / T_{s,\mathrm{fast}}``, onto the slow
baseband. The gain of the fast partition alone, from its inputs to that signal and evaluated
in this band, is a sufficient but conservative check: it ignores the attenuation the
continuous partition provides and rejects a boundary whose fast partition does not filter.
The check that accounts for the whole path assembles the composite of the continuous
partition, discretized at the fast rate, and the fast partition, and evaluates the gain from
all inputs of the composite to the boundary signal in the band. The held outputs of the slow
partition have to be among these inputs: the hold at the slow rate produces images of the
slow signal at multiples of its sampling frequency, which fall in this band and have to be
attenuated by the continuous partition.

```julia
dt1 = 0.1  # fast clock
dt2 = 0.3  # slow clock
k1 = ShiftIndex(Clock(dt1))
k2 = ShiftIndex(Clock(dt2))
@variables xc(t) d(t) w(t) z(t)
eqs = [
    D(xc) ~ -xc + Hold(z) + d                       # plant
    w(k1) ~ 0.5w(k1 - 1) + 0.5Sample(dt1)(xc)       # fast filter
    z(k2) ~ 0.9z(k2 - 1) - 0.1Sample(dt2)(Hold(w))  # slow controller
]
@named mr = System(eqs, t)
hlm = linearize_hybrid(mr, [d], [xc]; op = Dict(xc => 0.0, d => 0.0, w => 0.0, z => 0.0))
ic = hlm.continuous_index
ifast = findfirst(p -> p.Ts == dt1, hlm.partitions)
islow = findfirst(p -> p.Ts == dt2, hlm.partitions)
cont, fast = hlm.partitions[ic], hlm.partitions[ifast]
F = ss(fast.A, fast.B, fast.C, fast.D, fast.Ts)
Pf = c2d(ss(cont.A, cont.B, cont.C, cont.D), dt1)  # the plant at the fast rate
cont_to_fast = clock_boundary(hlm, ic, ifast)
fast_to_cont = clock_boundary(hlm, ifast, ic)      # no such signal in this model
fast_to_slow = clock_boundary(hlm, ifast, islow)
# The composite from all inputs of the plant (the disturbance and the held slow output) and
# the user-specified inputs of the fast partition to the signal entering the slow partition.
Gb = feedback(
    Pf, F; Y1 = cont_to_fast.outputs, U2 = cont_to_fast.inputs,
    Y2 = fast_to_cont.outputs, U1 = fast_to_cont.inputs,
    W1 = 1:Pf.nu, Z1 = Int[], W2 = 1:fast.nu_user, Z2 = fast_to_slow.outputs,
    pos_feedback = true
)
band = exp10.(range(log10(pi / dt2), log10(pi / dt1), length = 200))
below = exp10.(range(log10(pi / dt2) - 3, log10(pi / dt2), length = 200))
maximum(abs, freqresp(Gb, band)) / maximum(abs, freqresp(Gb, below))
```

For this model the ratio is about 0.06, whereas the fast filter alone attenuates the band by
less than a factor of two, so the attenuation of the signal entering the slow partition comes
from the plant. The content above the Nyquist frequency of the fast sampler is not
represented by the composite; it is checked on the continuous partition as above. With more
clocks, the signal entering a partition is checked on the composite of the continuous
partition and the partitions on faster clocks than the receiving one.

For multi-rate models, each discrete partition carries its own sample interval and the
connections span partitions on different clocks. Such loops can be analyzed either by
converting the discrete partitions to continuous time with `d2c` or by lifting them to a
common rate. In both cases the result is only meaningful for frequencies well below the
Nyquist frequency of every sampler in the loop, since the frequency content above the Nyquist
frequency of a sampler is folded by the sampling and is not represented by the individual
partition linearizations; the checks above indicate whether this content is small.

The Jacobians of each partition are computed from the generated functions of the compiled
partition with the AD backend `autodiff` (`AutoForwardDiff()` by default), so models calling
arbitrary Julia functions are supported. If differentiation of a partition fails, it is
retried with `fallback_autodiff` (`AutoFiniteDiff()` by default).

The operating point can be given as a dictionary or as a [`LinearizationOpPoint`](@ref)
wrapping a solution of the model and a time. It provides the values of the state variables of
every partition; the history variables of a discrete partition take the value of the variable
they are the history of, so the operating point is assumed to be stationary across ticks. The
value of a signal crossing a clock boundary is taken from the operating point if present and
otherwise evaluated from the operating point of the partition the signal originates from, so
that in the example above the held controller output follows from the plant state. The
initialization equations of the model that refer to the continuous partition only and
determine a value the operating point does not fix take part in the initialization of the
continuous partition; parameters bound to `missing` and constant states that such equations
determine therefore take their values from them. Remaining values that are not available
default to the guess of the symbol, or to zero, which is reported by a warning. When the
operating point is a `LinearizationOpPoint`, the values of the unknowns and parameters of
each partition are taken from the solution; the values the solution holds for observed
variables are not, since they would repeat the constraints the unknowns already satisfy.

Only periodic clocks are supported. Partitions on other clocks, as well as continuous and
discrete events, assertions and state machines of the model, are not accounted for and
produce a warning. Discrete partitions containing an algebraic loop within a single clock
cannot be linearized. Linearization is not defined with respect to Boolean-valued variables,
such as the state of a relay with hysteresis, and a partition containing them produces an
error naming the variables. An equation that is not differentiable at the operating point,
such as `sqrt(u)` at `u = 0`, gives non-finite entries in the matrices of its partition, which
is reported by a warning.

The analysis-point functions [`get_sensitivity`](@ref), [`get_comp_sensitivity`](@ref) and
[`get_looptransfer`](@ref) accept the keyword argument `hybrid = true` to return a
`HybridLinearization` instead of the continuous linearization, and `linearize_hybrid` accepts
analysis points as inputs and outputs in the same way as `linearize`.

## Docstrings

```@index
Pages = ["Linearization.md"]
```

```@docs; canonical = false
linearize
ModelingToolkit.linearize_symbolic
ModelingToolkit.linearization_function
ModelingToolkit.LinearizationProblem
ModelingToolkit.LinearizationOpPoint
linearize_hybrid
HybridLinearization
ClockPartitionLinearization
clock_boundary
```
