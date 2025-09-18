# Disturbance and input modeling modeling

Disturbances are often seen as external factors that influence a system. Modeling and simulation of such external influences is common in order to ensure that the plant and or control system can adequately handle or suppress these disturbances. Disturbance modeling is also integral to the problem of state estimation, indeed, modeling how disturbances affect the evolution of the state of the system is crucial in order to accurately estimate this state.

This tutorial will show how to model disturbances in ModelingToolkit as _disturbance inputs_. This involves demonstrating best practices that make it easy to use a single model to handle both disturbed and undisturbed systems, and making use of the model for both simulation and state estimation.

## A flexible component-based workflow

We will consider a simple system consisting of two inertias connected through a flexible shaft, such as a simple transmission system in a car. We start by modeling the plant _without any input signals_:

```@example DISTURBANCE_MODELING
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, Test
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Blocks
t = ModelingToolkit.t_nounits
D = ModelingToolkit.D_nounits

@mtkmodel SystemModel begin
    @parameters begin
        m1 = 1
        m2 = 1
        k = 10 # Spring stiffness
        c = 3  # Damping coefficient
    end
    @components begin
        inertia1 = Inertia(; J = m1, phi = 0, w = 0)
        inertia2 = Inertia(; J = m2, phi = 0, w = 0)
        spring = Spring(; c = k)
        damper = Damper(; d = c)
        torque = Torque(use_support = false)
    end
    @equations begin
        connect(torque.flange, inertia1.flange_a)
        connect(inertia1.flange_b, spring.flange_a, damper.flange_a)
        connect(inertia2.flange_a, spring.flange_b, damper.flange_b)
    end
end
```

Here, we have added a `torque` component that allows us to add a torque input to drive the system, but we have not connected any signal to it yet. We have not yet made any attempts at modeling disturbances, and this is deliberate, we will handle this later in order to make the plant model as generically useful as possible.

In order to simulate this system in the presence of disturbances, we must 1. Reason about how disturbances may affect the system, and 2. attach _disturbance inputs_ and _disturbance signals_ to the model. We distinguish between an _input_ and a _signal_ here, where we by _input_ mean an attachment point (connector) to which we may connect a _signal_, i.e., a time-varying function.

We create a new model that includes disturbance inputs and signals, and attach those to the already defined plant model. We assume that each of the two inertias can be affected by a disturbance torque, such as due to friction or an unknown load on the output inertia.

```@example DISTURBANCE_MODELING
@mtkmodel ModelWithInputs begin
    @components begin
        input_signal = Blocks.Sine(frequency = 1, amplitude = 1)
        disturbance_signal1 = Blocks.Step(height = -1, start_time = 2) # We add an input signal that equals zero by default so that it has no effect during normal simulation
        disturbance_signal2 = Blocks.Step(height = 2, start_time = 4)
        disturbance_torque1 = Torque(use_support = false)
        disturbance_torque2 = Torque(use_support = false)
        system_model = SystemModel()
    end
    @equations begin
        connect(input_signal.output, :u, system_model.torque.tau)
        connect(disturbance_signal1.output, :d1, disturbance_torque1.tau) # When we connect the input _signals_, we do so through an analysis point. This allows us to easily disconnect the input signals in situations when we do not need them. 
        connect(disturbance_signal2.output, :d2, disturbance_torque2.tau)
        connect(disturbance_torque1.flange, system_model.inertia1.flange_b)
        connect(disturbance_torque2.flange, system_model.inertia2.flange_b)
    end
end
```

This outer model, `ModelWithInputs`, contains two disturbance inputs, both of type `Torque`. It also contains three signal specifications, one for the control input and two for the corresponding disturbance inputs. Note how we added the disturbance torque inputs at this level of the model, but the control input was added inside the system model. This is a design choice that is up to the modeler, here, we consider the driving torque to be a fundamental part of the model that is always required to make use of it, while the disturbance inputs may be of interest only in certain situations, and we thus add them when needed. Since we have added not only input connectors, but also connected input signals to them, this model is complete and ready for simulation, i.e., there are no _unbound inputs_.

```@example DISTURBANCE_MODELING
@named model = ModelWithInputs() # Model with load disturbance
ssys = mtkcompile(model)
prob = ODEProblem(ssys, [], (0.0, 6.0))
sol = solve(prob, Tsit5())
using Plots
plot(sol)
```

A thing to note in the specification of `ModelWithInputs` is the presence of three [analysis points](https://docs.sciml.ai/ModelingToolkit/dev/tutorials/linear_analysis/#ModelingToolkit.AnalysisPoint), `:u`, `:d1`, and `:d2`. When signals are connected through an analysis point, we may at any time linearize the model as if the signals were not connected, i.e., as if the corresponding inputs were unbound. We may also use this to generate a julia function for the dynamics on the form ``f(x,u,p,t,w)`` where the input ``u`` and disturbance ``w`` may be provided as separate function arguments, as if the corresponding input signals were not present in the model. More details regarding this will be presented further below, here, we just demonstrate how we could linearize this system model from the inputs to the angular velocity of the inertias

```@example DISTURBANCE_MODELING
using ControlSystemsBase, ControlSystemsMTK # ControlSystemsMTK provides the high-level function named_ss and ControlSystemsBase provides the bodeplot function
P = named_ss(model, [ssys.u, ssys.d1, ssys.d2],
    [ssys.system_model.inertia1.w, ssys.system_model.inertia2.w])
bodeplot(P, plotphase = false)
```

It's worth noting at this point that the fact that we could connect disturbance outputs from outside of the plant-model definition was enabled by the fact that we used a component-based workflow, where the plant model had the appropriate connectors available. If the plant model had modeled the system using direct equations without connectors, this would not have been possible and the model would thus be significantly less flexible.

We summarize the findings so far as a number of best practices:

!!! tip "Best practices"
    
      - Use a component-based workflow to model the plant
      - If possible, model the plant without explicit disturbance inputs to make it as generic as possible
      - When disturbance inputs are needed, create a new model that includes the plant model and the disturbance inputs
      - Only add input _signals_ at the top level of the model, this applies to both control inputs and disturbance inputs.
      - Use analysis points to connect signals to inputs, this allows for easy disconnection of signals when needed, e.g., for linearization or function generation.

## Modeling for state estimation

In the example above, we constructed a model for _simulation_ of a disturbance affecting the system. When simulating, we connect an input signal of specified shape that simulates the disturbance, above, we used `Blocks.Step` as input signals. On the other hand, when performing state estimation, the exact shape of the disturbance is typically not known, we might only have some diffuse knowledge of the disturbance characteristics such as "varies smoothly", "makes sudden step changes" or "is approximately periodic with 24hr period". The encoding of such knowledge is commonly reasoned about in the frequency domain, where we specify a disturbance model as a dynamical system with a frequency response similar to the approximate spectrum of the disturbance. [For more details around this, see the in-depth tutorial notebook "How to tune a Kalman filter"](https://juliahub.com/pluto/editor.html?id=ad9ecbf9-bf83-45e7-bbe8-d2e5194f2240). Most algorithms for state estimation, such as a Kalman-filter like estimators, assume that disturbances are independent and identically distributed (i.i.d.). While seemingly restrictive at first glance, when combined with an appropriate disturbance models encoded as dynamical systems, this assumption still allows for a wide range of non i.i.d. disturbances to be modeled.

When modeling a system in MTK, we essentially (without considering algebraic equations for simplicity in exposition) construct a model of a dynamical system

```math
\begin{aligned}
\dot x &= f(x, p, t) \\
y &= g(x, p, t)
\end{aligned}
```

where ``x`` is the state, ``y`` are observed variables, ``p`` are parameters, and ``t`` is time. When using MTK, which variables constitute ``x`` and which are considered part of the output, ``y``, is up to the tool rather than the user, this choice is made by MTK during the call to `@mtkcompile` or the lower-level function `mtkcompile`.

If we further consider external inputs to the system, such as controlled input signals ``u`` and disturbance inputs ``w``, we can write the system as

```math
\begin{aligned}
\dot x &= f(x, u, p, t, w) \\
y &= g(x, u, p, t)
\end{aligned}
```

To make use of the model defined above for state estimation, we may want to generate a Julia function for the dynamics ``f`` and the output equations ``g`` that we can plug into, e.g., a nonlinear version of a Kalman filter or a particle filter, etc. MTK contains utilities to do this, namely, [`ModelingToolkit.generate_control_function`](@ref) and [`ModelingToolkit.build_explicit_observed_function`](@ref) (described in more details in ["Input output"](@ref inputoutput)). These functions take keyword arguments `disturbance_inputs` and `disturbance_argument`, that indicate which variables in the model are considered part of ``w``, and whether or not these variables are to be added as function arguments to ``f``, i.e., whether we have ``f(x, u, p, t)`` or ``f(x, u, p, t, w)``. If we do not include the disturbance inputs as function arguments, MTK will assume that the ``w`` variables are all zero, but any dynamics associated with these variables, such as disturbance models, will be included in the generated function. This allows a state estimator to estimate the state of the disturbance model, provided that this state is [observable](https://en.wikipedia.org/wiki/Observability) from the measured outputs of the system.

Below, we demonstrate

 1. How to add an integrating disturbance model
 2. how to generate the functions ``f`` and ``g`` for a typical nonlinear state estimator with explicit disturbance inputs

```@example DISTURBANCE_MODELING
@mtkmodel IntegratingDisturbance begin
    @variables begin
        x(t) = 0.0
        w(t) = 0.0, [disturbance = true, input = true]
    end
    @components begin
        input = RealInput()
        output = RealOutput()
    end
    @equations begin
        D(x) ~ w
        w ~ input.u
        output.u ~ x
    end
end

@mtkmodel SystemModelWithDisturbanceModel begin
    @components begin
        input_signal = Blocks.Sine(frequency = 1, amplitude = 1)
        disturbance_signal1 = Blocks.Constant(k = 0)
        disturbance_signal2 = Blocks.Constant(k = 0)
        disturbance_torque1 = Torque(use_support = false)
        disturbance_torque2 = Torque(use_support = false)
        disturbance_model = Blocks.Integrator()
        system_model = SystemModel()
    end
    @equations begin
        connect(input_signal.output, :u, system_model.torque.tau)
        connect(disturbance_signal1.output, :d1, disturbance_model.input)
        connect(disturbance_model.output, disturbance_torque1.tau)
        connect(disturbance_signal2.output, :d2, disturbance_torque2.tau)
        connect(disturbance_torque1.flange, system_model.inertia1.flange_b)
        connect(disturbance_torque2.flange, system_model.inertia2.flange_b)
    end
end

@named model_with_disturbance = SystemModelWithDisturbanceModel()
```

We demonstrate that this model is complete and can be simulated:

```@example DISTURBANCE_MODELING
ssys = mtkcompile(model_with_disturbance)
prob = ODEProblem(ssys, [], (0.0, 10.0))
sol = solve(prob, Tsit5())
using Test
@test SciMLBase.successful_retcode(sol)
```

but we may also generate the functions ``f`` and ``g`` for state estimation:

```@example DISTURBANCE_MODELING
inputs = [ssys.u]
disturbance_inputs = [ssys.d1, ssys.d2]
P = ssys.system_model
outputs = [P.inertia1.phi, P.inertia2.phi, P.inertia1.w, P.inertia2.w]

(f_oop, f_ip), x_sym,
p_sym,
io_sys = ModelingToolkit.generate_control_function(
    model_with_disturbance, inputs, disturbance_inputs; disturbance_argument = true)

g = ModelingToolkit.build_explicit_observed_function(
    io_sys, outputs; inputs)

op = ModelingToolkit.inputs(io_sys) .=> 0
x0 = ModelingToolkit.get_u0(io_sys, op)
p = MTKParameters(io_sys, op)
u = zeros(1) # Control input
w = zeros(length(disturbance_inputs)) # Disturbance input
@test f_oop(x0, u, p, t, w) == zeros(5)
@test g(x0, u, p, 0.0) == [0, 0, 0, 0]

# Non-zero disturbance inputs should result in non-zero state derivatives. We call `sort` since we do not generally know the order of the state variables
w = [1.0, 2.0]
@test sort(f_oop(x0, u, p, t, w)) == [0, 0, 0, 1, 2]
```

## Input signal library

The [`Blocks` module in ModelingToolkitStandardLibrary](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/blocks/) contains several predefined input signals, such as `Sine, Step, Ramp, Constant` etc., a few of which were used in the examples above. If you have an input signal represented as a sequence of samples, you may use an [`Interpolation` block](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/tutorials/input_component/), e.g., as `src = Interpolation(ConstantInterpolation, data, timepoints)`, see the docstring for a complete example.

## Disturbance-model library

There is no library explicitly constructed for disturbance modeling. Standard blocks from the [`Blocks` module in ModelingToolkitStandardLibrary](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/blocks/), such as `Integrator, TransferFunction, StateSpace`, can model any disturbance with rational spectrum. Examples of this includes disturbance models such as constants, piecewise constant, periodic, highpass, lowpass, and bandpass. For help with filter design, see [ControlSystems.jl: Filter-design](https://juliacontrol.github.io/ControlSystems.jl/stable/man/creating_systems/#Filter-design) and the interface package [ControlSystemsMTK.jl](https://juliacontrol.github.io/ControlSystemsMTK.jl/dev/). In the example above, we made use of `Blocks.Integrator`, which is a disturbance model suitable for disturbances dominated by low-frequency components, such as piecewise constant signals or slowly drifting signals.

## Further reading

To see full examples that perform state estimation with ModelingToolkit models, see the following resources:

  - [C codegen considered unnecessary: go directly to binary, do not pass C. Compilation of Julia code for deployment in model-based engineering](https://arxiv.org/abs/2502.01128)
  - [LowLevelParticleFiltersMTK.jl](https://github.com/baggepinnen/LowLevelParticleFiltersMTK.jl)

## Index

```@index
Pages = ["disturbance_modeling.md"]
```

```@autodocs; canonical = false
Modules = [ModelingToolkit]
Pages   = ["systems/analysis_points.jl"]
Order   = [:function, :type]
Private = false
```
