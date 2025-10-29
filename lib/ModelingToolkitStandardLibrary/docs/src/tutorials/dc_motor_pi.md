# DC Motor with PI-controller

In this example, a PI-controller is set up for speed control of a DC-motor. An equivalent circuit diagram is depicted below.

![DC-motor](https://user-images.githubusercontent.com/50108075/196108356-0e8605e3-61a9-4006-8559-786252e55928.png)

## Modeling and simulation

The electrical part consists of a resistance and inductance. The coupling between the electrical and rotational domain is done via an electro-motive force (EMF) component. The voltage across the EMF is proportional to the angular velocity and the current is proportional to the torque. On the mechanical side, viscous friction in, e.g., a bearing and the inertia of the shaft is modelled.

A PI-controller with anti-windup measure should be used as a speed controller. A simulation is performed to verify the tracking performance of the controller and the disturbance rejection capabilities.

First, the needed packages are imported and the parameters of the model defined.

```@example dc_motor_pi
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEq
using Plots
```

The actual model can now be composed.

```@example dc_motor_pi
@mtkmodel DCMotor begin
    @parameters begin
        R = 0.5, [description = "Armature resistance"] # Ohm
        L = 4.5e-3, [description = "Armature inductance"] # H
        k = 0.5, [description = "Motor constant"] # N.m/A
        J = 0.02, [description = "Inertia"] # kg.m²
        f = 0.01, [description = "Friction factor"] # N.m.s/rad
        tau_L_step = -0.3, [description = "Amplitude of the load torque step"] # N.m
    end
    @components begin
        ground = Ground()
        source = Voltage()
        ref = Blocks.Step(height = 1, start_time = 0)
        pi_controller = Blocks.LimPI(k = 1.1, T = 0.035, u_max = 10, Ta = 0.035)
        feedback = Blocks.Feedback()
        R1 = Resistor(R = R)
        L1 = Inductor(L = L)
        emf = EMF(k = k)
        fixed = Fixed()
        load = Torque()
        load_step = Blocks.Step(height = tau_L_step, start_time = 3)
        inertia = Inertia(J = J)
        friction = Damper(d = f)
        speed_sensor = SpeedSensor()
    end
    @equations begin
        connect(fixed.flange, emf.support, friction.flange_b)
        connect(emf.flange, friction.flange_a, inertia.flange_a)
        connect(inertia.flange_b, load.flange)
        connect(inertia.flange_b, speed_sensor.flange)
        connect(load_step.output, load.tau)
        connect(ref.output, feedback.input1)
        connect(speed_sensor.w, :y, feedback.input2)
        connect(feedback.output, pi_controller.err_input)
        connect(pi_controller.ctr_output, :u, source.V)
        connect(source.p, R1.p)
        connect(R1.n, L1.p)
        connect(L1.n, emf.p)
        connect(emf.n, source.n, ground.g)
    end
end

@named model = DCMotor()
nothing # hide
```

Now the model can be simulated. Typical rotational mechanical systems are described via `DAE`
(differential algebraic equations), however in this case, ModelingToolkit can simplify the model enough
so that it can be represented as a system of `ODEs` (ordinary differential equations).

```@example dc_motor_pi
sys = mtkcompile(model)
# Provide complete initial conditions for all state variables
u0 = Dict(
    sys.L1.i => 0.0,                  # Initial inductor current
    sys.inertia.w => 0.0,              # Initial angular velocity
    sys.inertia.phi => 0.0,            # Initial angle
    sys.pi_controller.int.x => 0.0    # Initial PI integrator state
)
prob = ODEProblem(sys, u0, (0, 6.0))
sol = solve(prob)

p1 = plot(sol.t, sol[sys.inertia.w], ylabel = "Angular Vel. in rad/s",
    label = "Measurement", title = "DC Motor with Speed Controller")
plot!(sol.t, sol[sys.ref.output.u], label = "Reference")
p2 = plot(sol.t, sol[sys.load.tau.u], ylabel = "Disturbance in Nm", label = "")
plot(p1, p2, layout = (2, 1))
```

## Closed-loop analysis

When implementing and tuning a control system in simulation, it is a good practice to analyze the closed-loop properties and verify robustness of the closed-loop with respect to, e.g., modeling errors. To facilitate this, we added two analysis points to the set of connections above, more specifically, we added the analysis points named `:y` and `:u` to the connections (for more details on analysis points, see [Linear Analysis](@ref))

```julia
connect(sys.speed_sensor.w, :y, sys.feedback.input2)
connect(sys.pi_controller.ctr_output, :u, sys.source.V)
```

one at the plant output (`:y`) and one at the plant input (`:u`). We may use these analysis points to calculate, e.g., sensitivity functions, illustrated below. Here, we calculate the sensitivity function $S(s)$ and the complimentary sensitivity function $T(s) = I - S(s)$, defined as

```math
\begin{aligned}
S(s) &= \dfrac{1}{I + P(s)C(s)} \\
T(s) &= \dfrac{P(s)C(s)}{I + P(s)C(s)}
\end{aligned}
```

```@example dc_motor_pi
using ControlSystemsBase
# Get sensitivity function
matrices_S,
simplified_sys_S = Blocks.get_sensitivity(
    model, :y, op = Dict(unknowns(sys) .=> 0.0))
So = ss(matrices_S...) |> minreal # The output-sensitivity function as a StateSpace system
# Get complementary sensitivity function
matrices_T,
simplified_sys_T = Blocks.get_comp_sensitivity(
    model, :y, op = Dict(unknowns(sys) .=> 0.0))
To = ss(matrices_T...)# The output complementary sensitivity function as a StateSpace system
bodeplot([So, To], label = ["S" "T"], plot_title = "Sensitivity functions",
    plotphase = false)
```

Similarly, we may compute the loop-transfer function and plot its Nyquist curve

```@example dc_motor_pi
matrices_L,
simplified_sys_L = Blocks.get_looptransfer(
    model, :y, op = Dict(unknowns(sys) .=> 0.0))
L = -ss(matrices_L...) # The loop-transfer function as a StateSpace system. The negative sign is to negate the built-in negative feedback
Ms, ωMs = hinfnorm(So) # Compute the peak of the sensitivity function to draw a circle in the Nyquist plot
nyquistplot(L, label = "\$L(s)\$", ylims = (-2.5, 0.5), xlims = (-1.2, 0.1),
    Ms_circles = Ms)
```
