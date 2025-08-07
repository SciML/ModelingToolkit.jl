#=
This file implements and tests a typical workflow for state estimation with disturbance models
The primary subject of the tests is the analysis-point features and the
analysis-point specific method for `generate_control_function`.
=#
using ModelingToolkit, OrdinaryDiffEqTsit5, LinearAlgebra, Test
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Blocks
import NonlinearSolve
using ModelingToolkit: connect
# using Plots

using ModelingToolkit: t_nounits as t, D_nounits as D

indexof(sym, syms) = findfirst(isequal(sym), syms)

## Build the system model ======================================================
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

# The addition of disturbance inputs relies on the fact that the plant model has been constructed using connectors, we use these to connect the disturbance inputs from outside the plant-model definition
@mtkmodel ModelWithInputs begin
    @components begin
        input_signal = Blocks.Sine(frequency = 1, amplitude = 1)
        disturbance_signal1 = Blocks.Constant(k = 0) # We add an input signal that equals zero by default so that it has no effect during normal simulation
        disturbance_signal2 = Blocks.Constant(k = 0)
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

@named model = ModelWithInputs() # Model with load disturbance
ssys = mtkcompile(model)
prob = ODEProblem(ssys, [], (0.0, 10.0))
sol = solve(prob, Tsit5())
# plot(sol)

##
using ControlSystemsBase, ControlSystemsMTK
cmodel = complete(model)
P = cmodel.system_model
lsys = named_ss(
    model, [:u, :d1], [P.inertia1.phi, P.inertia2.phi, P.inertia1.w, P.inertia2.w])

##
# If we now want to add a disturbance model, we cannot do that since we have already connected a constant to the disturbance input, we thus create a new wrapper model with inputs

s = tf("s")
dist(; name) = System(1 / s; name)

@mtkmodel SystemModelWithDisturbanceModel begin
    @components begin
        input_signal = Blocks.Sine(frequency = 1, amplitude = 1)
        disturbance_signal1 = Blocks.Constant(k = 0)
        disturbance_signal2 = Blocks.Constant(k = 0)
        disturbance_torque1 = Torque(use_support = false)
        disturbance_torque2 = Torque(use_support = false)
        disturbance_model = dist()
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
# ssys = mtkcompile(open_loop(model_with_disturbance, :d)) # Open loop worked, but it's a bit awkward that we have to use it here
# lsys2 = named_ss(model_with_disturbance, [:u, :d1],
# [P.inertia1.phi, P.inertia2.phi, P.inertia1.w, P.inertia2.w])
ssys = mtkcompile(model_with_disturbance)
prob = ODEProblem(ssys, [], (0.0, 10.0))
sol = solve(prob, Tsit5())
@test SciMLBase.successful_retcode(sol)
# plot(sol)

## 
# Now we only have an integrating disturbance affecting inertia1, what if we want both integrating and direct Gaussian? We'd need a "PI controller" disturbancemodel. If we add the disturbance model (s+1)/s we get the integrating and non-integrating noises being correlated which is fine, it reduces the dimensions of the sigma point by 1.

dist3(; name) = System(ss(1 + 10 / s, balance = false); name)

@mtkmodel SystemModelWithDisturbanceModel begin
    @components begin
        input_signal = Blocks.Sine(frequency = 1, amplitude = 1)
        disturbance_signal1 = Blocks.Constant(k = 0)
        disturbance_signal2 = Blocks.Constant(k = 0)
        disturbance_torque1 = Torque(use_support = false)
        disturbance_torque2 = Torque(use_support = false)
        disturbance_model = dist3()
        system_model = SystemModel()

        y = Blocks.Add()
        angle_sensor = AngleSensor()
        output_disturbance = Blocks.Constant(k = 0)
    end
    @equations begin
        connect(input_signal.output, :u, system_model.torque.tau)
        connect(disturbance_signal1.output, :d1, disturbance_model.input)
        connect(disturbance_model.output, disturbance_torque1.tau)
        connect(disturbance_signal2.output, :d2, disturbance_torque2.tau)
        connect(disturbance_torque1.flange, system_model.inertia1.flange_b)
        connect(disturbance_torque2.flange, system_model.inertia2.flange_b)

        connect(system_model.inertia1.flange_b, angle_sensor.flange)
        connect(angle_sensor.phi, y.input1)
        connect(output_disturbance.output, :dy, y.input2)
    end
end

@named model_with_disturbance = SystemModelWithDisturbanceModel()
# ssys = mtkcompile(open_loop(model_with_disturbance, :d)) # Open loop worked, but it's a bit awkward that we have to use it here
# lsys3 = named_ss(model_with_disturbance, [:u, :d1],
#     [P.inertia1.phi, P.inertia2.phi, P.inertia1.w, P.inertia2.w])
ssys = mtkcompile(model_with_disturbance)
prob = ODEProblem(ssys, [], (0.0, 10.0))
sol = solve(prob, Tsit5())
@test SciMLBase.successful_retcode(sol)
# plot(sol)

## Generate function for an augmented Unscented Kalman Filter =====================
# temp = open_loop(model_with_disturbance, :d)
outputs = [P.inertia1.phi, P.inertia2.phi, P.inertia1.w, P.inertia2.w]
f, x_sym,
p_sym,
io_sys = ModelingToolkit.generate_control_function(
    model_with_disturbance, [:u], [:d1, :d2, :dy], split = false)

f, x_sym,
p_sym,
io_sys = ModelingToolkit.generate_control_function(
    model_with_disturbance, [:u], [:d1, :d2, :dy],
    disturbance_argument = true, split = false)

measurement = ModelingToolkit.build_explicit_observed_function(
    io_sys, outputs, inputs = ModelingToolkit.inputs(io_sys)[1:1])
measurement2 = ModelingToolkit.build_explicit_observed_function(
    io_sys, [io_sys.y.output.u], inputs = ModelingToolkit.inputs(io_sys)[1:1],
    disturbance_inputs = ModelingToolkit.inputs(io_sys)[2:end],
    disturbance_argument = true)

op = ModelingToolkit.inputs(io_sys) .=> 0
x0 = ModelingToolkit.get_u0(io_sys, op)
p = ModelingToolkit.get_p(io_sys, op)
x = zeros(5)
u = zeros(1)
d = zeros(3)
@test f[1](x, u, p, t, d) == zeros(5)
@test measurement(x, u, p, 0.0) == [0, 0, 0, 0]
@test measurement2(x, u, p, 0.0, d) == [0]

# Add to the integrating disturbance input
d = [1, 0, 0]
@test sort(f[1](x, u, p, 0.0, d)) == [0, 0, 0, 1, 1] # Affects disturbance state and one velocity
@test measurement2(x, u, p, 0.0, d) == [0]

d = [0, 1, 0]
@test sort(f[1](x, u, p, 0.0, d)) == [0, 0, 0, 0, 1] # Affects one velocity
@test measurement(x, u, p, 0.0) == [0, 0, 0, 0]
@test measurement2(x, u, p, 0.0, d) == [0]

d = [0, 0, 1]
@test sort(f[1](x, u, p, 0.0, d)) == [0, 0, 0, 0, 0] # Affects nothing
@test measurement(x, u, p, 0.0) == [0, 0, 0, 0]
@test measurement2(x, u, p, 0.0, d) == [1] # We have now disturbed the output

## Further downstream tests that the functions generated above actually have the properties required to use for state estimation
# 
# using LowLevelParticleFilters, SeeToDee
# Ts = 0.001
# discrete_dynamics = SeeToDee.Rk4(f_oop2, Ts)
# nx = length(x_sym)
# nu = 1
# nw = 2
# ny = length(outputs)
# R1 = Diagonal([1e-5, 1e-5])
# R2 = 0.1 * I(ny)
# op = ModelingToolkit.inputs(io_sys) .=> 0
# x0, p = ModelingToolkit.get_u0_p(io_sys, op, op)
# d0 = LowLevelParticleFilters.SimpleMvNormal(x0, 10.0I(nx))
# measurement_model = UKFMeasurementModel{Float64, false, false}(measurement, R2; nx, ny)
# kf = UnscentedKalmanFilter{false, false, true, false}(
#     discrete_dynamics, measurement_model, R1, d0; nu, Ts, p)

# tvec = 0:Ts:sol.t[end]
# u = vcat.(Array(sol(tvec, idxs = P.torque.tau.u)))
# y = collect.(eachcol(Array(sol(tvec, idxs = outputs)) .+ 1e-2 .* randn.()))

# inds = 1:5805
# res = forward_trajectory(kf, u, y)

# plot(res, size = (1000, 1000));
# plot!(sol, idxs = x_sym, sp = (1:nx)', l = :dash);
