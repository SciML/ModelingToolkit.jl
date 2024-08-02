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
@variables x(t)=0 y(t)=0 u(t)=0 r(t)=0
@parameters kp = 1

eqs = [u ~ kp * (r - y) # P controller
       D(x) ~ -x + u    # First-order plant
       y ~ x]           # Output equation

@named sys = ODESystem(eqs, t) # Do not call @mtkbuild when linearizing
matrices, simplified_sys = linearize(sys, [r], [y]) # Linearize from r to y
matrices
```

The named tuple `matrices` contains the matrices of the linear statespace representation, while `simplified_sys` is an `ODESystem` that, among other things, indicates the unknown variable order in the linear system through

```@example LINEARIZE
using ModelingToolkit: inputs, outputs
[unknowns(simplified_sys); inputs(simplified_sys); outputs(simplified_sys)]
```

!!! note "Inputs must be unconnected"
    
    The model above has 4 variables but only three equations, there is no equation specifying the value of `r` since `r` is an input. This means that only unbalanced models can be linearized, or in other words, models that are balanced and can be simulated _cannot_ be linearized. To learn more about this, see [How to linearize a ModelingToolkit model (YouTube)](https://www.youtube.com/watch?v=-XOux-2XDGI&t=395s). Also see [ModelingToolkitStandardLibrary: Linear analysis](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/linear_analysis/) for utilities that make linearization of completed models easier.

!!! note "Un-simplified system"
    
    Linearization expects `sys` to be un-simplified, i.e., `structural_simplify` or `@mtkbuild` should not be called on the system before linearizing.

## Operating point

The operating point to linearize around can be specified with the keyword argument `op` like this: `op = Dict(x => 1, r => 2)`. The operating point may include specification of unknown variables, input variables and parameters. For variables that are not specified in `op`, the default value specified in the model will be used if available, if no value is specified, an error is thrown.

## Batch linearization and algebraic variables

If linearization is to be performed around multiple operating points, the simplification of the system has to be carried out a single time only. To facilitate this, the lower-level function [`ModelingToolkit.linearization_function`](@ref) is available. This function further allows you to obtain separate Jacobians for the differential and algebraic parts of the model. For ODE models without algebraic equations, the statespace representation above is available from the output of `linearization_function` as `A, B, C, D = f_x, f_u, h_x, h_u`.

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

[ModelingToolkitStandardLibrary](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/) contains a set of [tools for more advanced linear analysis](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/linear_analysis/). These can be used to make it easier to work with and analyze causal models, such as control and signal-processing systems.

Also see [ControlSystemsMTK.jl](https://juliacontrol.github.io/ControlSystemsMTK.jl/dev/) for an interface to [ControlSystems.jl](https://github.com/JuliaControl/ControlSystems.jl) that contains tools for linear analysis and frequency-domain analysis.

## Docstrings

```@index
Pages = ["Linearization.md"]
```

```@docs
linearize
ModelingToolkit.linearize_symbolic
ModelingToolkit.linearization_function
```
