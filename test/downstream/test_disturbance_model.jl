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

@testset "Symbolic known disturbances with an analysis point" begin
    @variables x(t) = 0.0 source(t) [output = true] control(t) [input = true] w(t) [input = true]
    @named sys = System([source ~ 0, connect(source, :u, control), D(x) ~ -x + control + w], t)

    @testset "Control input $(typeof(input))" for input in (:u, sys.u)
        (; f, io_sys) = ModelingToolkit.generate_control_function(
            sys, input; known_disturbance_inputs = [w], split = false
        )
        op = ModelingToolkit.inputs(io_sys) .=> 0.0
        p = ModelingToolkit.get_p(io_sys, op)
        @test f[1]([0.0], [0.0], p, 0.0, [2.0]) == [2.0]
    end
end

@testset "Marked known disturbance at an analysis point" begin
    @variables x(t) control(t) [input = true] control_source(t) [output = true]
    @variables w(t) [input = true, disturbance = true] disturbance_source(t) [output = true]
    @variables unknown_disturbance(t) [disturbance = true]
    @named sys = System(
        [
            control_source ~ 0, disturbance_source ~ 0,
            connect(control_source, :u, control), connect(disturbance_source, :d, w),
            D(x) ~ -x + 5control - 3w + unknown_disturbance,
        ], t
    )

    (; f, io_sys) = ModelingToolkit.generate_control_function(
        sys, sys.u; known_disturbance_inputs = [sys.d], split = false
    )
    p = ModelingToolkit.get_p(io_sys, ModelingToolkit.inputs(io_sys) .=> 0.0)
    @test f[1]([2.0], [3.0], p, 0.0, [4.0]) == [1.0]
end

@testset "Array-valued known disturbance at an analysis point (split = $split)" for split in (true, false)
    @variables x(t) control(t) [input = true] control_source(t) [output = true]
    @variables (w(t))[1:2] [input = true] (disturbance_source(t))[1:2] [output = true]
    @named sys = System(
        [
            control_source ~ 0, disturbance_source ~ [0.0, 0.0],
            connect(control_source, :u, control), connect(disturbance_source, :d, w),
            D(x) ~ -x + 5control - 3w[1] + 7w[2],
        ], t
    )

    @testset "$interface" for (interface, args, kwargs) in (
            ("keyword", (sys.u,), (; known_disturbance_inputs = [sys.d])),
            ("positional", (sys.u, sys.d), (; disturbance_argument = true)),
        )
        (; f, ps, io_sys) = ModelingToolkit.generate_control_function(
            sys, args...; split, kwargs...
        )
        @test isempty(ps)
        p = ModelingToolkit.get_p(io_sys, ModelingToolkit.inputs(io_sys) .=> 0.0)
        @test f[1]([2.0], [3.0], p, 0.0, [4.0, 5.0]) == [36.0]
        @test f[1]([2.0], [3.0], p, 0.0, [5.0, 5.0]) == [33.0]
        @test f[1]([2.0], [3.0], p, 0.0, [4.0, 6.0]) == [43.0]
    end
end

energy(v) = sum(abs2, v)
@register_symbolic energy(v::AbstractVector)

@testset "Array-valued unknown disturbance at an analysis point (split = $split)" for split in (true, false)
    @variables x(t) c(t) [input = true] control_source(t) [output = true]
    @variables (w(t))[1:2] [input = true] (load_source(t))[1:2] [output = true]
    @named sys = System(
        [
            control_source ~ 0, load_source ~ [0.0, 0.0],
            connect(control_source, :control_ap, c), connect(load_source, :load_ap, w),
            D(x) ~ 5c + energy(w),
        ], t
    )

    (; f, io_sys) = ModelingToolkit.generate_control_function(sys, sys.control_ap, sys.load_ap; split)
    p = ModelingToolkit.get_p(io_sys, ModelingToolkit.inputs(io_sys) .=> 0.0)
    @test f[1]([0.0], [3.0], p, 0.0) == [15.0]
end

@testset "Mixed disturbances at array elements (split = $split)" for split in (true, false)
    @variables x(t) c(t) [input = true] control_source(t) [output = true]
    @variables (w(t))[1:2] [input = true, disturbance = true] load_source(t) [output = true]
    @testset "Known elements $indices" for (indices, load, indexed_expected, array_expected) in (
            ((1,), [4.0], 17.0, 31.0),
            ((2,), [4.0], -23.0, 31.0),
            ((2, 1), [4.0, 6.0], -5.0, 67.0),
        )
        @testset "$consumer" for (consumer, rhs, control, expected) in (
                ("Indexed", 5c + 3w[1] - 7w[2], [1.0], indexed_expected),
                ("Whole array", 5c + energy(w), [3.0], array_expected),
            )
            @named sys = System(
                [
                    control_source ~ 0, load_source ~ 0,
                    connect(control_source, :control_ap, c), connect(load_source, :load_ap, w[first(indices)]),
                    D(x) ~ rhs,
                ], t
            )

            known_disturbances = [sys.load_ap, (w[i] for i in Base.tail(indices))...]
            (; f, io_sys) = ModelingToolkit.generate_control_function(
                sys, sys.control_ap; known_disturbance_inputs = known_disturbances, split
            )
            p = ModelingToolkit.get_p(io_sys, ModelingToolkit.inputs(io_sys) .=> 0.0)
            @test f[1]([0.0], control, p, 0.0, load) == [expected]
            dx = [NaN]
            f[2](dx, [0.0], control, p, 0.0, load)
            @test dx == [expected]
        end
    end
end

indexof(sym, syms) = findfirst(isequal(sym), syms)

## Build the system model ======================================================
@component function SystemModel(; name, m1 = 1, m2 = 1, k = 10, c = 3)
    pars = @parameters begin
        m1 = m1
        m2 = m2
        k = k # Spring stiffness
        c = c  # Damping coefficient
    end

    systems = @named begin
        inertia1 = Inertia(; J = m1, phi = 0, w = 0)
        inertia2 = Inertia(; J = m2, phi = 0, w = 0)
        spring = Spring(; c = k)
        damper = Damper(; d = c)
        torque = Torque(use_support = false)
    end

    vars = @variables begin
    end

    equations = Equation[
        connect(torque.flange, inertia1.flange_a),
        connect(inertia1.flange_b, spring.flange_a, damper.flange_a),
        connect(inertia2.flange_a, spring.flange_b, damper.flange_b),
    ]

    return System(equations, t, vars, pars; name, systems)
end

# The addition of disturbance inputs relies on the fact that the plant model has been constructed using connectors, we use these to connect the disturbance inputs from outside the plant-model definition
@component function ModelWithInputs(; name)
    pars = @parameters begin
    end

    systems = @named begin
        input_signal = Blocks.Sine(frequency = 1, amplitude = 1)
        disturbance_signal1 = Blocks.Constant(k = 0) # We add an input signal that equals zero by default so that it has no effect during normal simulation
        disturbance_signal2 = Blocks.Constant(k = 0)
        disturbance_torque1 = Torque(use_support = false)
        disturbance_torque2 = Torque(use_support = false)
        system_model = SystemModel()
    end

    vars = @variables begin
    end

    equations = Equation[
        connect(input_signal.output, :u, system_model.torque.tau),
        connect(disturbance_signal1.output, :d1, disturbance_torque1.tau), # When we connect the input _signals_, we do so through an analysis point. This allows us to easily disconnect the input signals in situations when we do not need them.
        connect(disturbance_signal2.output, :d2, disturbance_torque2.tau),
        connect(disturbance_torque1.flange, system_model.inertia1.flange_b),
        connect(disturbance_torque2.flange, system_model.inertia2.flange_b),
    ]

    return System(equations, t, vars, pars; name, systems)
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
    model, [:u, :d1], [P.inertia1.phi, P.inertia2.phi, P.inertia1.w, P.inertia2.w]
)

##
# If we now want to add a disturbance model, we cannot do that since we have already connected a constant to the disturbance input, we thus create a new wrapper model with inputs

s = tf("s")
dist(; name) = System(1 / s; name)

@component function SystemModelWithDisturbanceModel(; name)
    pars = @parameters begin
    end

    systems = @named begin
        input_signal = Blocks.Sine(frequency = 1, amplitude = 1)
        disturbance_signal1 = Blocks.Constant(k = 0)
        disturbance_signal2 = Blocks.Constant(k = 0)
        disturbance_torque1 = Torque(use_support = false)
        disturbance_torque2 = Torque(use_support = false)
        disturbance_model = dist()
        system_model = SystemModel()
    end

    vars = @variables begin
    end

    equations = Equation[
        connect(input_signal.output, :u, system_model.torque.tau),
        connect(disturbance_signal1.output, :d1, disturbance_model.input),
        connect(disturbance_model.output, disturbance_torque1.tau),
        connect(disturbance_signal2.output, :d2, disturbance_torque2.tau),
        connect(disturbance_torque1.flange, system_model.inertia1.flange_b),
        connect(disturbance_torque2.flange, system_model.inertia2.flange_b),
    ]

    return System(equations, t, vars, pars; name, systems)
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

@component function SystemModelWithDisturbanceModel(; name)
    pars = @parameters begin
    end

    systems = @named begin
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

    vars = @variables begin
    end

    equations = Equation[
        connect(input_signal.output, :u, system_model.torque.tau),
        connect(disturbance_signal1.output, :d1, disturbance_model.input),
        connect(disturbance_model.output, disturbance_torque1.tau),
        connect(disturbance_signal2.output, :d2, disturbance_torque2.tau),
        connect(disturbance_torque1.flange, system_model.inertia1.flange_b),
        connect(disturbance_torque2.flange, system_model.inertia2.flange_b),

        connect(system_model.inertia1.flange_b, angle_sensor.flange),
        connect(angle_sensor.phi, y.input1),
        connect(output_disturbance.output, :dy, y.input2),
    ]

    return System(equations, t, vars, pars; name, systems)
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
    model_with_disturbance, [:u], [:d1, :d2, :dy], split = false
)

f, x_sym,
    p_sym,
    io_sys = ModelingToolkit.generate_control_function(
    model_with_disturbance, [:u], [:d1, :d2, :dy],
    disturbance_argument = true, split = false
)

@testset "Known analysis-point disturbances ($(eltype(disturbance_inputs)))" for (control_inputs, disturbance_inputs) in (
        ([:u], [:d1, :d2]),
        ([model_with_disturbance.u], [model_with_disturbance.d1, model_with_disturbance.d2]),
    )
    (; f, dvs, io_sys) = ModelingToolkit.generate_control_function(
        model_with_disturbance, control_inputs;
        known_disturbance_inputs = disturbance_inputs, split = false
    )
    @test length(ModelingToolkit.inputs(io_sys)) == length(control_inputs) + length(disturbance_inputs)
    op = ModelingToolkit.inputs(io_sys) .=> 0
    x0 = ModelingToolkit.get_u0(io_sys, op)
    p = ModelingToolkit.get_p(io_sys, op)
    u = zeros(length(control_inputs))
    @test f[1](x0, u, p, 0.0, zeros(length(disturbance_inputs))) == zeros(length(x0))
    inertia1_velocity = io_sys.system_model.inertia1.w
    inertia2_velocity = io_sys.system_model.inertia2.w
    integrator_state = io_sys.disturbance_model.x[1]
    expected_d1 = [isequal(x, inertia1_velocity) || isequal(x, integrator_state) for x in dvs]
    expected_d2 = [isequal(x, inertia2_velocity) for x in dvs]
    @test f[1](x0, u, p, 0.0, [1.0, 0.0]) == expected_d1
    @test f[1](x0, u, p, 0.0, [0.0, 2.0]) == 2 .* expected_d2
    @test f[1](x0, u, p, 0.0, [1.0, 2.0]) == expected_d1 + 2 .* expected_d2
end

measurement = ModelingToolkit.build_explicit_observed_function(
    io_sys, outputs, inputs = ModelingToolkit.inputs(io_sys)[1:1]
)
measurement2 = ModelingToolkit.build_explicit_observed_function(
    io_sys, [io_sys.y.output.u], inputs = ModelingToolkit.inputs(io_sys)[1:1],
    disturbance_inputs = ModelingToolkit.inputs(io_sys)[2:end],
    disturbance_argument = true
)
# Test new interface with known_disturbance_inputs (equivalent to above)
measurement3 = ModelingToolkit.build_explicit_observed_function(
    io_sys, [io_sys.y.output.u], inputs = ModelingToolkit.inputs(io_sys)[1:1],
    known_disturbance_inputs = ModelingToolkit.inputs(io_sys)[2:end]
)

op = ModelingToolkit.inputs(io_sys) .=> 0
x0 = ModelingToolkit.get_u0(io_sys, op)
p = ModelingToolkit.get_p(io_sys, op)
x = zeros(5)
u = zeros(1)
d = zeros(3)
@test f[1](x, u, p, t, d) == zeros(5)
@test measurement(x, u, p, 0.0) == [0, 0, 0, 0]
@test measurement2(x, u, p, 0.0, d) == [0]
@test measurement3(x, u, p, 0.0, d) == [0]  # Test new interface

# Add to the integrating disturbance input
d = [1, 0, 0]
@test sort(f[1](x, u, p, 0.0, d)) == [0, 0, 0, 1, 1] # Affects disturbance state and one velocity
@test measurement2(x, u, p, 0.0, d) == [0]
@test measurement3(x, u, p, 0.0, d) == [0]  # Test new interface

d = [0, 1, 0]
@test sort(f[1](x, u, p, 0.0, d)) == [0, 0, 0, 0, 1] # Affects one velocity
@test measurement(x, u, p, 0.0) == [0, 0, 0, 0]
@test measurement2(x, u, p, 0.0, d) == [0]
@test measurement3(x, u, p, 0.0, d) == [0]  # Test new interface

d = [0, 0, 1]
@test sort(f[1](x, u, p, 0.0, d)) == [0, 0, 0, 0, 0] # Affects nothing
@test measurement(x, u, p, 0.0) == [0, 0, 0, 0]
@test measurement2(x, u, p, 0.0, d) == [1] # We have now disturbed the output
@test measurement3(x, u, p, 0.0, d) == [1]  # Test new interface

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
