# Linear Analysis

Linear analysis refers to the process of linearizing a nonlinear model and analysing the resulting linear dynamical system. To facilitate linear analysis, ModelingToolkit provides the concept of an [`AnalysisPoint`](@ref), which can be inserted in-between two causal blocks (such as those from `ModelingToolkitStandardLibrary.Blocks` sub module). Once a model containing analysis points is built, several operations are available:

  - [`get_sensitivity`](@ref) get the [sensitivity function (wiki)](https://en.wikipedia.org/wiki/Sensitivity_(control_systems)), $S(s)$, as defined in the field of control theory.
  - [`get_comp_sensitivity`](@ref) get the complementary sensitivity function $T(s) : S(s)+T(s)=1$.
  - [`get_looptransfer`](@ref) get the (open) loop-transfer function where the loop starts and ends in the analysis point. For a typical simple feedback connection with a plant $P(s)$ and a controller $C(s)$, the loop-transfer function at the plant output is $P(s)C(s)$.
  - [`linearize`](@ref) can be called with two analysis points denoting the input and output of the linearized system.
  - [`open_loop`](@ref) return a new (nonlinear) system where the loop has been broken in the analysis point, i.e., the connection the analysis point usually implies has been removed.

An analysis point can be created explicitly using the constructor [`AnalysisPoint`](@ref), or automatically when connecting two causal components using `connect`:

```julia
connect(comp1.output, :analysis_point_name, comp2.input)
```

A single output can also be connected to multiple inputs:

```julia
connect(comp1.output, :analysis_point_name, comp2.input, comp3.input, comp4.input)
```

!!! warning "Causality"
    
    Analysis points are *causal*, i.e., they imply a directionality for the flow of information. The order of the connections in the connect statement is thus important, i.e., `connect(out, :name, in)` is different from `connect(in, :name, out)`.

The directionality of an analysis point can be thought of as an arrow in a block diagram, where the name of the analysis point applies to the arrow itself.

```
┌─────┐         ┌─────┐
│     │  name   │     │
│  out├────────►│in   │
│     │         │     │
└─────┘         └─────┘
```

This is signified by the name being the middle argument to `connect`.

Of the above mentioned functions, all except for [`open_loop`](@ref) return the output of [`ModelingToolkit.linearize`](@ref), which is

```julia
matrices, simplified_sys = linearize(...)
# matrices = (; A, B, C, D)
```

i.e., `matrices` is a named tuple containing the matrices of a linear state-space system on the form

```math
\begin{aligned}
\dot x &= Ax + Bu\\
y &= Cx + Du
\end{aligned}
```

## Example

The following example builds a simple closed-loop system with a plant $P$ and a controller $C$. Two analysis points are inserted, one before and one after $P$. We then derive a number of sensitivity functions and show the corresponding code using the package ControlSystemBase.jl

```@example LINEAR_ANALYSIS
using ModelingToolkitStandardLibrary.Blocks, ModelingToolkit
@named P = FirstOrder(k = 1, T = 1) # A first-order system with pole in -1
@named C = Gain(-1)             # A P controller
t = ModelingToolkit.get_iv(P)
eqs = [connect(P.output, :plant_output, C.input)  # Connect with an automatically created analysis point called :plant_output
       connect(C.output, :plant_input, P.input)]
sys = ODESystem(eqs, t, systems = [P, C], name = :feedback_system)

matrices_S = get_sensitivity(sys, :plant_input)[1] # Compute the matrices of a state-space representation of the (input)sensitivity function.
matrices_T = get_comp_sensitivity(sys, :plant_input)[1]
```

Continued linear analysis and design can be performed using ControlSystemsBase.jl.
We create `ControlSystemsBase.StateSpace` objects using

```@example LINEAR_ANALYSIS
using ControlSystemsBase, Plots
S = ss(matrices_S...)
T = ss(matrices_T...)
bodeplot([S, T], lab = ["S" "" "T" ""], plot_title = "Bode plot of sensitivity functions",
    margin = 5Plots.mm)
```

The sensitivity functions obtained this way should be equivalent to the ones obtained with the code below

```@example LINEAR_ANALYSIS_CS
using ControlSystemsBase
P = tf(1.0, [1, 1]) |> ss
C = 1                      # Negative feedback assumed in ControlSystems
S = sensitivity(P, C)      # or feedback(1, P*C)
T = comp_sensitivity(P, C) # or feedback(P*C)
```

We may also derive the loop-transfer function $L(s) = P(s)C(s)$ using

```@example LINEAR_ANALYSIS
matrices_L = get_looptransfer(sys, :plant_output)[1]
L = ss(matrices_L...)
```

which is equivalent to the following with ControlSystems

```@example LINEAR_ANALYSIS_CS
L = P * (-C) # Add the minus sign to build the negative feedback into the controller
```

To obtain the transfer function between two analysis points, we call `linearize`

```@example LINEAR_ANALYSIS
using ModelingToolkit # hide
matrices_PS = linearize(sys, :plant_input, :plant_output)[1]
```

this particular transfer function should be equivalent to the linear system `P(s)S(s)`, i.e., equivalent to

```@example LINEAR_ANALYSIS_CS
feedback(P, C)
```

### Obtaining transfer functions

A statespace system from [ControlSystemsBase](https://juliacontrol.github.io/ControlSystems.jl/stable/man/creating_systems/) can be converted to a transfer function using the function `tf`:

```@example LINEAR_ANALYSIS_CS
tf(S)
```

## Gain and phase margins

Further linear analysis can be performed using the [analysis methods from ControlSystemsBase](https://juliacontrol.github.io/ControlSystems.jl/stable/lib/analysis/). For example, calculating the gain and phase margins of a system can be done using

```@example LINEAR_ANALYSIS_CS
margin(P)
```

(they are infinite for this system). A Nyquist plot can be produced using

```@example LINEAR_ANALYSIS_CS
nyquistplot(P)
```

## Index

```@index
Pages = ["linear_analysis.md"]
```

```@autodocs
Modules = [ModelingToolkit]
Pages   = ["systems/analysis_points.jl"]
Order   = [:function, :type]
Private = false
```
