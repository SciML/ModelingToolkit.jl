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

The `linearize` function expects the user to specify the inputs ``u`` and the outputs ``u`` using the syntax shown in the example below:

## Example
```@example LINEARIZE
using ModelingToolkit
@variables t x(t)=0 y(t)=0 u(t)=0 r(t)=0
@parameters kp = 1
D = Differential(t)

eqs = [u ~ kp * (r - y) # P controller
       D(x) ~ -x + u    # First-order plant
       y ~ x]           # Output equation

@named sys = ODESystem(eqs, t)
matrices, simplified_sys = linearize(sys, [r], [y]) # Linearize from r to y
matrices
```
The named tuple `matrices` contains the matrices of the linear statespace representation, while `simplified_sys` is an `ODESystem` that, amongst other things, indicates the state order in the linear system through
```@example LINEARIZE
using ModelingToolkit: inputs, outputs
[states(simplified_sys); inputs(simplified_sys); outputs(simplified_sys)]
```

## Operating point
The operating point to linearize around can be specified with the keyword argument `op` like this: `op = Dict(x => 1, r => 2)`. 

## Batch linearization and algebraic variables
If linearization is to be performed around multiple operating points, the simplification of the system has to be carried out a single time only. To facilitate this, the lower-level function [`ModelingToolkit.linearization_function`](@ref) is available. This function further allows you to obtain separate Jacobians for the differential and algebraic parts of the model. For ODE models without algebraic equations, the statespace representation above is available from the output of `linearization_function` as `A, B, C, D = f_x, f_u, h_x, h_u`.


## Input derivatives
Physical systems are always *proper*, i.e., they do not differentiate causal inputs. However, ModelingToolkit allows you to model non-proper systems, such as inverse models, and may sometimes fail to find a realization of a proper system on proper form. In these situations, `linearize` may throw an error mentioning
```
Input derivatives appeared in expressions (-g_z\g_u != 0)
```
This means that to simulate this system, some order of derivatives of the input is required. To allow `linearize` to proceed in this situation, one may pass the keyword argument `allow_input_derivatives = true`, in which case the resulting model will have twice as many inputs, ``2n_u``, where the last ``n_u`` inputs correspond to ``\dot u``. 

If the modeled system is actually proper (but MTK failed to find a proper realization), further numerical simplification can be applied to the resulting statespace system to obtain a proper form. Such simplification is currently available in the experimental package [ControlSystemsMTK](https://github.com/baggepinnen/ControlSystemsMTK.jl#internals-transformation-of-non-proper-models-to-proper-statespace-form).


## Tools for linear analysis
[ModelingToolkitStandardLibrary](http://mtkstdlib.sciml.ai/dev/API/linear_analysis/) contains a set of [tools for more advanced linear analysis](http://mtkstdlib.sciml.ai/dev/API/linear_analysis/). These can be used to make it easier to work with and analyze causal models, such as control and signal-processing systems. 

```@index
Pages = ["Linearization.md"]
```

```@docs
linearize
ModelingToolkit.linearization_function
```