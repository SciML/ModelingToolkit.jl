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

@named duffing = System(eqs, t, systems = [y, u], defaults = [u.u => 0])

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

## Docstrings

```@index
Pages = ["Linearization.md"]
```

```@docs; canonical = false
linearize
ModelingToolkit.linearize_symbolic
ModelingToolkit.linearization_function
ModelingToolkit.LinearizationProblem
```
