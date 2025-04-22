# [Input output](@id inputoutput)

An input-output system is a system on the form

```math
\begin{aligned}
M \dot x &= f(x, u, p, t) \\
y &= g(x, u, p, t)
\end{aligned}
```

where ``x`` is the state, ``u`` is the input and ``y`` is an output (in some contexts called an _observed variables_ in MTK).

While many uses of ModelingToolkit for simulation do not require the user to think about inputs and outputs (IO), there are certain situations in which handling IO explicitly may be important, such as

  - Linearization
  - Control design
  - System identification
  - FMU export
  - Real-time simulation with external data inputs
  - Custom interfacing with other simulation tools

This documentation page lists utilities that are useful for working with inputs and outputs in ModelingToolkit.

## Generating a dynamics function with inputs, ``f``

ModelingToolkit can generate the dynamics of a system, the function ``M\dot x = f(x, u, p, t)`` above, such that the user can pass not only the state ``x`` and parameters ``p`` but also an external input ``u``. To this end, the function [`ModelingToolkit.generate_control_function`](@ref) exists.

This function takes a vector of variables that are to be considered inputs, i.e., part of the vector ``u``. Alongside returning the function ``f``, [`ModelingToolkit.generate_control_function`](@ref) also returns the chosen state realization of the system after simplification. This vector specifies the order of the state variables ``x``, while the user-specified vector `u` specifies the order of the input variables ``u``.

!!! note "Un-simplified system"
    
    This function expects `sys` to be un-simplified, i.e., `structural_simplify` or `@mtkbuild` should not be called on the system before passing it into this function. `generate_control_function` calls a special version of `structural_simplify` internally.

### Example:

The following example implements a simple first-order system with an input `u` and state `x`. The function `f` is generated using `generate_control_function`, and the function `f` is then tested with random input and state values.

```@example inputoutput
using ModelingToolkit
import ModelingToolkit: t_nounits as t, D_nounits as D
@variables x(t)=0 u(t)=0 y(t)
@parameters k = 1
eqs = [D(x) ~ -k * (x + u)
       y ~ x]

@named sys = ODESystem(eqs, t)
f, x_sym, ps = ModelingToolkit.generate_control_function(sys, [u], simplify = true);
nothing # hide
```

We can inspect the state realization chosen by MTK

```@example inputoutput
x_sym
```

as expected, `x` is chosen as the state variable.

```@example inputoutput
using Test # hide
@test isequal(x_sym[], x) # hide
@test isequal(ps, [k]) # hide
nothing  # hide
```

Now we can test the generated function `f` with random input and state values

```@example inputoutput
p = [1]
x = [rand()]
u = [rand()]
@test f[1](x, u, p, 1) ≈ -p[] * (x + u) # Test that the function computes what we expect D(x) = -k*(x + u)
```

## Generating an output function, ``g``

ModelingToolkit can also generate a function that computes a specified output of a system, the function ``y = g(x, u, p, t)`` above. This is done using the function [`build_explicit_observed_function`](@ref). When generating an output function, the user must specify the output variable(s) of interest, as well as any inputs if inputs are relevant to compute the output.

The order of the user-specified output variables determines the order of the output vector ``y``.

## Input-output variable metadata

See [Symbolic Metadata](@ref symbolic_metadata). Metadata specified when creating variables is not directly used by any of the functions above, but the user can use the accessor functions `ModelingToolkit.inputs(sys)` and `ModelingToolkit.outputs(sys)` to obtain all variables with such metadata for passing to the functions above. The presence of this metadata is not required for any IO functionality and may be omitted.

## Linearization

See [Linearization](@ref linearization).

## Docstrings

```@index
Pages = ["InputOutput.md"]
```

```@docs
ModelingToolkit.generate_control_function
ModelingToolkit.build_explicit_observed_function
```
