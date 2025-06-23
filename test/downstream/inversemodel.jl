using ModelingToolkit
using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using SymbolicIndexingInterface
using Test
using ControlSystemsMTK: tf, ss, get_named_sensitivity, get_named_comp_sensitivity
using ModelingToolkit: t_nounits as t, D_nounits as D
# ==============================================================================
## Mixing tank
# This tests a common workflow in control engineering, the use of an inverse-based
# feedforward model. Such a model differentiates "inputs", exercising the dummy-derivative functionality of ModelingToolkit. We also test linearization and computation of sensitivity functions
# for such models.
# ==============================================================================

connect = ModelingToolkit.connect;
rc = 0.25 # Reference concentration 

@mtkmodel MixingTank begin
    @parameters begin
        c0 = 0.8, [description = "Nominal concentration"]
        T0 = 308.5, [description = "Nominal temperature"]
        a1 = 0.2674
        a21 = 1.815
        a22 = 0.4682
        b = 1.5476
        k0 = 1.05e14
        ϵ = 34.2894
    end
    @variables begin
        gamma(t), [description = "Reaction speed"]
        xc(t) = c0, [description = "Concentration"]
        xT(t) = T0, [description = "Temperature"]
        xT_c(t), [description = "Cooling temperature"]
    end
    @components begin
        T_c = RealInput()
        c = RealOutput()
        T = RealOutput()
    end
    begin
        τ0 = 60
        wk0 = k0 / c0
        wϵ = ϵ * T0
        wa11 = a1 / τ0
        wa12 = c0 / τ0
        wa13 = c0 * a1 / τ0
        wa21 = a21 / τ0
        wa22 = a22 * T0 / τ0
        wa23 = T0 * (a21 - b) / τ0
        wb = b / τ0
    end
    @equations begin
        gamma ~ xc * wk0 * exp(-wϵ / xT)
        D(xc) ~ -wa11 * xc - wa12 * gamma + wa13
        D(xT) ~ -wa21 * xT + wa22 * gamma + wa23 + wb * xT_c
        xc ~ c.u
        xT ~ T.u
        xT_c ~ T_c.u
    end
end
begin
    Ftf = tf(1, [(100), 1])^2
    Fss = ss(Ftf)
    # Create an MTK-compatible constructor 
    function RefFilter(; name)
        sys = System(Fss; name)
        "Compute initial state that yields y0 as output"
        empty!(ModelingToolkit.get_defaults(sys))
        return sys
    end
end
@mtkmodel InverseControlledTank begin
    begin
        c0 = 0.8    #  "Nominal concentration
        T0 = 308.5 #  "Nominal temperature
        x10 = 0.42
        x20 = 0.01
        u0 = -0.0224
        c_start = c0 * (1 - x10) # Initial concentration
        T_start = T0 * (1 + x20) # Initial temperature
        c_high_start = c0 * (1 - 0.72) # Reference concentration
        T_c_start = T0 * (1 + u0) # Initial cooling temperature
    end
    @components begin
        ref = Constant(k = 0.25) # Concentration reference
        ff_gain = Gain(k = 1) # To allow turning ff off
        controller = PI(gainPI.k = 10, T = 500)
        tank = MixingTank(xc = c_start, xT = T_start, c0 = c0, T0 = T0)
        inverse_tank = MixingTank(xc = nothing, xT = T_start, c0 = c0, T0 = T0)
        feedback = Feedback()
        add = Add()
        filter = RefFilter()
        noise_filter = FirstOrder(k = 1, T = 1, x = T_start)
        # limiter = Gain(k=1)
        limiter = Limiter(y_max = 370, y_min = 250) # Saturate the control input 
    end
    @equations begin
        connect(ref.output, :r, filter.input)
        connect(filter.output, inverse_tank.c)
        connect(inverse_tank.T_c, ff_gain.input)
        connect(ff_gain.output, :uff, limiter.input)
        connect(limiter.output, add.input1)
        connect(controller.ctr_output, :u, add.input2)
        connect(add.output, :u_tot, tank.T_c)
        connect(inverse_tank.T, feedback.input1)
        connect(tank.T, :y, noise_filter.input)
        connect(noise_filter.output, feedback.input2)
        connect(feedback.output, :e, controller.err_input)
    end
end;
@named model = InverseControlledTank()
ssys = mtkcompile(model)
cm = complete(model)

op = Dict(
    cm.filter.y.u => 0.8 * (1 - 0.42),
    D(cm.filter.y.u) => 0
)
tspan = (0.0, 1000.0)
# https://github.com/SciML/ModelingToolkit.jl/issues/2786
prob = ODEProblem(ssys, op, tspan)
sol = solve(prob, Rodas5P())

@test SciMLBase.successful_retcode(sol)

# plot(sol, idxs=[model.tank.xc, model.tank.xT, model.controller.ctr_output.u], layout=3, sp=[1 2 3])
# hline!([prob[cm.ref.k]], label="ref", sp=1)

@test sol(tspan[2], idxs = cm.tank.xc)≈getp(prob, cm.ref.k)(prob) atol=1e-2 # Test that the inverse model led to the correct reference

# we need to provide `op` so the initialization system knows what to hold constant
# the values don't matter
Sf, simplified_sys = get_sensitivity_function(model, :y; op); # This should work without providing an operating opint containing a dummy derivative
x = state_values(Sf)
p = parameter_values(Sf)
# If this somehow passes, mention it on
# https://github.com/SciML/ModelingToolkit.jl/issues/2786
matrices1 = Sf(x, p, 0)
matrices2, _ = get_sensitivity(model, :y; op); # Test that we get the same result when calling the higher-level API
@test matrices1.f_x ≈ matrices2.A[1:6, 1:6]
nsys = get_named_sensitivity(model, :y; op) # Test that we get the same result when calling an even higher-level API
@test matrices2.A ≈ nsys.A

# Test the same thing for comp sensitivities

# This should work without providing an operating opint containing a dummy derivative
Sf, simplified_sys = get_comp_sensitivity_function(model, :y; op);
x = state_values(Sf)
p = parameter_values(Sf)
# If this somehow passes, mention it on
# https://github.com/SciML/ModelingToolkit.jl/issues/2786
matrices1 = Sf(x, p, 0)
# Test that we get the same result when calling the higher-level API
matrices2, _ = get_comp_sensitivity(model, :y; op)
@test matrices1.f_x ≈ matrices2.A[1:6, 1:6]
# Test that we get the same result when calling an even higher-level API
nsys = get_named_comp_sensitivity(model, :y; op)
@test matrices2.A ≈ nsys.A

@testset "Issue #3319" begin
    op1 = Dict(
        cm.filter.y.u => 0.8 * (1 - 0.42),
        cm.tank.xc => 0.65
    )

    op2 = Dict(
        cm.filter.y.u => 0.8 * (1 - 0.42),
        cm.tank.xc => 0.45
    )

    output = :y
    # we need to provide `op` so the initialization system knows which
    # values to hold constant
    lin_fun, ssys = get_sensitivity_function(model, output; op = op1)
    matrices1, extras1 = linearize(ssys, lin_fun, op = op1)
    matrices2, extras2 = linearize(ssys, lin_fun, op = op2)
    @test extras1.x != extras2.x
    S1f = ss(matrices1...)
    S2f = ss(matrices2...)
    @test S1f != S2f

    matrices1, ssys = get_sensitivity(model, output; op = op1)
    matrices2, ssys = get_sensitivity(model, output; op = op2)
    S1 = ss(matrices1...)
    S2 = ss(matrices2...)
    @test S1 != S2

    @test S1 == S1f
    @test S2 == S2f
end
