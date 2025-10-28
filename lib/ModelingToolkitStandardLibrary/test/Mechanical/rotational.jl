using ModelingToolkitStandardLibrary.Mechanical.Rotational,
      ModelingToolkit, OrdinaryDiffEq,
      Test
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEq: ReturnCode.Success

# using Plots

@testset "two inertias" begin
    @mtkmodel TwoInertia begin
        @components begin
            fixed = Fixed()
            inertia1 = Inertia(J = 2) # this one is fixed
            spring = Spring(c = 1e4)
            damper = Damper(d = 10)
            inertia2 = Inertia(J = 2, phi = pi / 2)
        end
        @equations begin
            connect(fixed.flange, inertia1.flange_b)
            connect(inertia1.flange_b, spring.flange_a, damper.flange_a)
            connect(spring.flange_b, damper.flange_b, inertia2.flange_a)
        end
    end

    @mtkcompile sys = TwoInertia()

    prob = ODEProblem(sys, [D(D(sys.inertia2.phi)) => 0], (0, 10.0))
    sol1 = solve(prob, Rodas4())
    @test SciMLBase.successful_retcode(sol1)

    prob = DAEProblem(
        sys, D.(unknowns(sys)) .=> prob.f(sol1.u[1], prob.p, 0.0), (0, 10.0))
    dae_sol = solve(prob, DFBDF())
    @test SciMLBase.successful_retcode(dae_sol)
    @test all(dae_sol[sys.inertia1.w] .== 0)
    @test dae_sol[sys.inertia2.w][end]≈0 atol=1e-3 # all energy has dissipated

    @mtkmodel WithSpringDamper begin
        @components begin
            fixed = Fixed()
            inertia1 = Inertia(J = 2) # this one is fixed
            springdamper = SpringDamper(; c = 1e4, d = 10)
            inertia2 = Inertia(J = 2, phi = pi / 2)
        end
        @equations begin
            connect(fixed.flange, inertia1.flange_b)
            connect(inertia1.flange_b, springdamper.flange_a)
            connect(springdamper.flange_b, inertia2.flange_a)
        end
    end

    @mtkcompile sys = WithSpringDamper()

    prob = ODEProblem(sys, [D(D(sys.inertia2.phi)) => 0], (0, 10.0))
    sol2 = solve(prob, Rodas4())
    @test SciMLBase.successful_retcode(sol2)
    @test sol2(0:1:10, idxs = sys.inertia2.w).u≈sol1(0:1:10, idxs = sys.inertia2.w).u atol=1e-3

    # Plots.plot(sol; vars=[inertia1.w, inertia2.w])
end

@testset "two inertias with driving torque" begin
    @mtkmodel TwoInertiasWithDrivingTorque begin
        @structural_parameters begin
            amplitude = 10 # Amplitude of driving torque
            frequency = 5 # Frequency of driving torque
            J_motor = 0.1 # Motor inertia
        end

        @components begin
            fixed = Fixed()
            torque = Torque(; use_support = true)
            inertia1 = Inertia(J = 2, phi = pi / 2)
            spring = Rotational.Spring(c = 1e4)
            damper = Damper(d = 10)
            inertia2 = Inertia(J = 4)
            sine = Blocks.Sine(amplitude = amplitude, frequency = frequency)
        end

        @equations begin
            connect(sine.output, torque.tau)
            connect(torque.support, fixed.flange)
            connect(torque.flange, inertia1.flange_a)
            connect(inertia1.flange_b, spring.flange_a, damper.flange_a)
            connect(spring.flange_b, damper.flange_b, inertia2.flange_a)
        end
    end

    @mtkcompile sys = TwoInertiasWithDrivingTorque()
    deqs = [eq.lhs => eq.rhs for eq in equations(sys)]
    prob = DAEProblem(
        sys, [deqs;
              D(D(sys.inertia2.phi)) => 1.0; sys.spring.flange_b.phi => 0.0], (0, 10.0))
    sol = solve(prob, DFBDF())
    @test SciMLBase.successful_retcode(sol)

    prob = ODEProblem(
        sys, [D(D(sys.inertia2.phi)) => 0.0, sys.spring.flange_b.phi => 0.0], (0, 1.0))
    sol = solve(prob, Rodas4())
    @test SciMLBase.successful_retcode(sol)

    # exact opposite oscillation with smaller amplitude J2 = 2*J1 and with an offset.
    @test all(isapprox.(
        sol[sys.inertia1.w], -sol[sys.inertia2.w] * 2 .+ sol[sys.inertia1.w][1], atol = 1))
    @test all(sol[sys.torque.flange.tau] .== -sol[sys.sine.output.u]) # torque source is equal to negative sine

    ## Test with constant torque source
    @mtkmodel TwoInertiasWitConstantTorque begin
        @components begin
            fixed = Fixed()
            torque = ConstantTorque(use_support = true, tau_constant = 1)
            inertia1 = Inertia(J = 2, phi = pi / 2)
            spring = Rotational.Spring(c = 1e4)
            damper = Damper(d = 10)
            inertia2 = Inertia(J = 4)
        end
        @equations begin
            connect(torque.support, fixed.flange)
            connect(torque.flange, inertia1.flange_a)
            connect(inertia1.flange_b, spring.flange_a, damper.flange_a)
            connect(spring.flange_b, damper.flange_b, inertia2.flange_a)
        end
    end

    @mtkcompile sys = TwoInertiasWitConstantTorque()

    prob = ODEProblem(
        sys, [D(D(sys.inertia2.phi)) => 1.0, sys.spring.flange_b.phi => 0.0], (0, 10.0))
    sol = solve(prob, Rodas4())
    @test SciMLBase.successful_retcode(sol)
    @test sol(sol.t[end], idxs = sys.inertia1.w)≈sol(sol.t[end], idxs = sys.inertia2.w) rtol=0.1 # both inertias have same angular velocity after initial transient
end

# see: https://doc.modelica.org/Modelica%204.0.0/Resources/helpWSM/Modelica/Modelica.Mechanics.Rotational.Examples.First.html
@testset "first example" begin
    @mtkmodel FirstExample begin
        @structural_parameters begin
            amplitude = 10 # Amplitude of driving torque
            frequency = 5 # Frequency of driving torque
            J_motor = 0.1 # Motor inertia
            J_load = 2 # Load inertia
            ratio = 10 # Gear ratio
            damping = 10 # Damping in bearing of gear
        end

        @components begin
            fixed = Fixed()
            torque = Torque(use_support = true)
            inertia1 = Inertia(J = J_motor)
            idealGear = IdealGear(ratio = ratio, use_support = true)
            inertia2 = Inertia(J = 2)
            spring = Spring(c = 1e4)
            inertia3 = Inertia(J = J_load)
            damper = Damper(d = damping)
            sine = Blocks.Sine(amplitude = amplitude, frequency = frequency)
        end

        @equations begin
            connect(inertia1.flange_b, idealGear.flange_a)
            connect(idealGear.flange_b, inertia2.flange_a)
            connect(inertia2.flange_b, spring.flange_a)
            connect(spring.flange_b, inertia3.flange_a)
            connect(damper.flange_a, inertia2.flange_b)
            connect(damper.flange_b, fixed.flange)
            connect(sine.output, torque.tau)
            connect(torque.support, fixed.flange)
            connect(idealGear.support, fixed.flange)
            connect(torque.flange, inertia1.flange_a)
        end
    end

    @mtkcompile sys = FirstExample()
    prob = ODEProblem(
        sys, [sys.inertia3.w => 0.0, sys.spring.flange_a.phi => 0.0], (0, 1.0))
    sol = solve(prob, Rodas4())
    @test SciMLBase.successful_retcode(sol)
    # Plots.plot(sol; vars=[inertia2.w, inertia3.w])
end

@testset "Stick-Slip" begin
    @mtkmodel VelocityProfile begin
        @components begin
            sine = Blocks.Sine(amplitude = 10, frequency = 0.1)
            dz = Blocks.DeadZone(u_max = 2)
            lim = Blocks.Limiter(y_max = 6)
            output = Blocks.RealOutput()
        end
        @equations begin
            connect(sine.output, dz.input)
            connect(dz.output, lim.input)
            connect(lim.output, output)
        end
    end

    @mtkmodel StickSlip begin
        @components begin
            fixed = Fixed()
            spring = Spring(c = 6.5)
            damper = Damper(d = 0.01)
            inertia = Inertia(J = 0.0001)
            friction = RotationalFriction(f = 0.001, tau_c = 20, w_brk = 0.06035,
                tau_brk = 25)
            vel_profile = VelocityProfile()
            source = Speed()
            angle_sensor = AngleSensor()
        end

        @equations begin
            connect(vel_profile.output, source.w_ref)
            connect(source.flange, friction.flange_a)
            connect(friction.flange_b, inertia.flange_a)
            connect(inertia.flange_b, spring.flange_a, damper.flange_a)
            connect(spring.flange_b, damper.flange_b, fixed.flange)
            connect(angle_sensor.flange, inertia.flange_a)
        end
    end

    @mtkcompile sys = StickSlip()
    prob = DAEProblem(sys,
        [D.(unknowns(sys)) .=> 0.0;
         [sys.inertia.flange_b.tau => 0.0; unknowns(sys) .=> 0.0...]],
        (0, 10.0))

    sol = solve(prob, DFBDF())
    @test SciMLBase.successful_retcode(sol)
    @test sol[sys.angle_sensor.phi.u] == sol[sys.inertia.flange_a.phi]

    # p1 = Plots.plot(sol; vars=[inertia.flange_a.phi, source.phi], title="Angular Position", labels=["Inertia" "Source"], ylabel="Angle in rad")
    # p2 = Plots.plot(sol; vars=[friction.w_rel], title="Rel. Angular Velocity of Friction", label="", ylabel="Angular Velocity in rad/s")
    # Plots.plot(p1, p2, layout=(2, 1))
    # Plots.savefig("stick_slip.png")

    # Plots.scatter(sol[friction.w], sol[friction.tau], label="")
end

@testset "sensors" begin
    @mtkmodel Sensors begin
        @components begin
            fixed = Fixed()
            inertia1 = Inertia(J = 2) # this one is fixed
            spring = Spring(c = 1e4)
            damper = Damper(d = 10)
            inertia2 = Inertia(J = 2, phi = pi / 2)
            speed_sensor = SpeedSensor()
            torque_sensor = TorqueSensor()
            rel_speed_sensor = RelSpeedSensor()
        end

        @equations begin
            connect(fixed.flange, inertia1.flange_b, rel_speed_sensor.flange_b)
            connect(inertia1.flange_b, torque_sensor.flange_a)
            connect(torque_sensor.flange_b, spring.flange_a, damper.flange_a,
                speed_sensor.flange, rel_speed_sensor.flange_a)
            connect(spring.flange_b, damper.flange_b, inertia2.flange_a)
        end
    end

    @mtkcompile sys = Sensors()

    prob = ODEProblem(sys, [D(D(sys.inertia2.phi)) => 0.0], (0, 10.0))
    sol = solve(prob, Rodas4())
    @test SciMLBase.successful_retcode(sol)
    @test all(sol[sys.inertia1.w] .== 0)
    @test all(sol[sys.inertia1.w] .== sol[sys.speed_sensor.w.u])
    @test sol[sys.inertia2.w][end]≈0 atol=1e-3 # all energy has dissipated
    @test all(sol[sys.rel_speed_sensor.w_rel.u] .== sol[sys.speed_sensor.w.u])
    @test all(sol[sys.torque_sensor.tau.u] .== -sol[sys.inertia1.flange_b.tau])

    prob = DAEProblem(
        sys, D.(unknowns(sys)) .=> prob.f(sol.u[1], prob.p, 0.0), (0, 10.0))
    sol = solve(prob, DFBDF())
    @test SciMLBase.successful_retcode(sol)
    @test all(sol[sys.inertia1.w] .== 0)
    @test all(sol[sys.inertia1.w] .== sol[sys.speed_sensor.w.u])
    @test sol[sys.inertia2.w][end]≈0 atol=1e-3 # all energy has dissipated
    @test all(sol[sys.rel_speed_sensor.w_rel.u] .== sol[sys.speed_sensor.w.u])
    @test all(sol[sys.torque_sensor.tau.u] .== -sol[sys.inertia1.flange_b.tau])

    # Plots.plot(sol; vars=[inertia1.w, inertia2.w])
end

@testset "Position" begin
    @mtkmodel TestPosition begin
        @components begin
            pos = Rotational.Position(exact = true, f_crit = 500)
            input = Blocks.Sine(frequency = 1, amplitude = 1)
            inertia = Rotational.Inertia(J = 1)
        end
        @equations begin
            connect(input.output, pos.phi_ref)
            connect(pos.flange, inertia.flange_a)
        end
    end
    @mtkcompile sys = TestPosition()
    prob = ODEProblem(sys, [], (0, 10.0))
    sol = solve(prob, Tsit5())
    @test SciMLBase.successful_retcode(sol)
    tv = 0:0.01:10
    @test sol(tv, idxs = sys.inertia.phi).u≈sin.(2pi .* tv) atol=1e-12

    @mtkmodel TestPosition begin
        @components begin
            pos = Rotational.Position(exact = false, f_crit = 500)
            input = Blocks.Sine(frequency = 1, amplitude = 1)
            inertia = Rotational.Inertia(J = 1)
        end
        @equations begin
            connect(input.output, pos.phi_ref)
            connect(pos.flange, inertia.flange_a)
        end
    end
    @mtkcompile sys = TestPosition()
    prob = ODEProblem(sys, [
            sys.inertia.phi => 0,
            sys.inertia.w => 0
        ], (0, 10.0))
    sol = solve(prob, Tsit5())
    @test SciMLBase.successful_retcode(sol)
    tv = 0:0.01:10
    @test sol(tv, idxs = sys.inertia.phi).u≈sin.(2pi .* tv) atol=1e-1
end
