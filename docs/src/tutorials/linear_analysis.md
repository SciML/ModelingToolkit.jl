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

An analysis point of a built model is referred to by property access on the system in
which it was created, `sys.analysis_point_name`, and the functions above accept the
[`AnalysisPoint`](@ref) this returns. A point belonging to a subsystem is reached through
the hierarchy, `sys.inner.analysis_point_name`. The system must be the unsimplified and
uncompleted model, since both `mtkcompile` and `complete` remove the analysis points.

!!! warning "Old style"

    Referring to an analysis point by its `Symbol` name, as in
    `get_sensitivity(sys, :plant_input)`, is the old style and is retained only for
    backwards compatibility. It will be deprecated.

Of the above mentioned functions, all except for [`open_loop`](@ref) return the output of [`ModelingToolkit.linearize`](@ref), which is

```julia
matrices, simplified_sys = linearize(_...)
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
using ModelingToolkit: t_nounits as t

@named P = FirstOrder(k = 1, T = 1) # A first-order system with pole in -1
@named C = Gain(-1)             # A P controller

eqs = [connect(P.output, :plant_output, C.input)  # Connect with an automatically created analysis point called :plant_output
       connect(C.output, :plant_input, P.input)]
sys = System(eqs, t, systems = [P, C], name = :feedback_system)

matrices_S = get_sensitivity(sys, sys.plant_input)[1] # Compute the matrices of a state-space representation of the (input)sensitivity function.
matrices_T = get_comp_sensitivity(sys, sys.plant_input)[1]
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
matrices_L = get_looptransfer(sys, sys.plant_output)[1]
L = ss(matrices_L...)
```

which is equivalent to the following with ControlSystems

```@example LINEAR_ANALYSIS_CS
L = P * (-C) # Add the minus sign to build the negative feedback into the controller
```

To obtain the transfer function between two analysis points, we call `linearize`

```@example LINEAR_ANALYSIS
using ModelingToolkit # hide
matrices_PS = linearize(sys, sys.plant_input, sys.plant_output)[1]
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

## Operating Point of Signals Disconnected by Loop Openings
When a connection is broken via `loop_openings`, the input variables that were driven by the connection are disconnected. Their values are part of the operating point at which the linearization is computed, which is significant for systems that are nonlinear in these variables.

Semantics: A loop opening changes the input-output structure of the linearization, not the operating point. A disconnected variable becomes a parameter of the system, not an additional input of the linearization, so it does not appear as an extra column of the $B$ or $D$ matrices. Its value is determined by, in order of precedence,

 1. its entry in the operating point: a value, e.g., `op = Dict(u => value)`; the value `missing`, which makes the initialization solve for the variable from additional conditions in the operating point, for example a steady-state condition `D(x) => 0` for a differential state variable `x`; or an expression,
 2. its initial condition in the system, e.g., `@variables u(t) = value`, unless the linearization is performed around a solution with [`LinearizationOpPoint`](@ref), which ignores the initial conditions of the system,
 3. otherwise, the value of the signal that drove it before the connection was broken, at the operating point of the loop-closed system. The initialization includes the equation `u ~ driver` in place of the broken connection, so it solves the same problem as the initialization of the loop-closed system. With a [`LinearizationOpPoint`](@ref), this is the value the signal has in the solution at each time point.

A guess value of the variable does not determine its value. The same rules apply to the input that replaces the broken connection in [`get_looptransfer`](@ref), and to an input analysis point that is also listed in `loop_openings`. The value of a disconnected signal is fixed when the linearization function is built, so it must be given in the operating point passed to [`linearization_function`](@ref) or to the analysis function that builds it; `linearize(sys, lin_fun; op)` rejects values for signals that the linearization function holds at the value of their driver.

Motivation:

  - Loop openings are mostly used to extract a part of a feedback system at the operating point of the closed loop, e.g., the plant, the controller, or a loop-transfer function for the computation of stability margins. The linearizations of the parts combine into the linearization of the closed loop only if all parts are linearized at the same operating point. Holding the disconnected signals at their loop-closed values ensures this without any entries in the operating point, both at a single operating point and along a solution.
  - Setting disconnected signals to zero changes the operating point of the components that are nonlinear in them, and can yield linearizations at points the system never attains.
  - A loop opening can also be used to disconnect a signal in order to linearize at a different operating point, for example to remove a reference step that saturates a controller at the initial time. The new value of the signal is then a choice of the user, and it is given in the operating point or as an initial condition of the disconnected variable.

## Index

```@index
Pages = ["linear_analysis.md"]
```

```@autodocs; canonical = false
Modules = [ModelingToolkit]
Pages   = ["systems/analysis_points.jl"]
Order   = [:function, :type]
Private = false
```
