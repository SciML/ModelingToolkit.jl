using ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkit: t_nounits as t, D_nounits as D

using ModelingToolkitStandardLibrary.Blocks
import ModelingToolkitStandardLibrary.Mechanical.Translational as TV
import ModelingToolkitStandardLibrary.Mechanical.TranslationalPosition as TP

@testset "Free" begin
    function TestSystem(; name)
        systems = @named begin
            acc = TV.Acceleration(false)
            a = Constant(; k = -10)
            mass = TV.Mass(; m = 100)
            free = TV.Free()
        end

        eqs = [connect(a.output, acc.a)
               connect(mass.flange, acc.flange, free.flange)]

        System(eqs, t, [], []; name, systems)
    end

    @named system = TestSystem()
    s = complete(system)
    sys = mtkcompile(system)
    prob = ODEProblem(sys, [s.mass.s => 0], (0, 0.1))
    sol = solve(prob, Rosenbrock23())

    @test sol[s.mass.flange.v][end]≈-0.1 * 10 atol=1e-3
    @test sol[s.free.f][end] ≈ 100 * 10
end

@testset "Spring, Damper, Mass, Fixed" begin
    @named dv = TV.Damper(d = 1)
    @named dp = TP.Damper(d = 1)

    @named sv = TV.Spring(k = 1)
    @named sp = TP.Spring(k = 1, l = 1)

    @named bv = TV.Mass(m = 1)
    @named bp = TP.Mass(m = 1, v = 1, s = 3)

    @named gv = TV.Fixed()
    @named gp = TP.Fixed(s_0 = 1)

    function simplify_and_solve(
            damping, spring, body, ground; initialization_eqs = Equation[])
        eqs = [connect(spring.flange_a, body.flange, damping.flange_a)
               connect(spring.flange_b, damping.flange_b, ground.flange)]

        @named model = System(eqs, t; systems = [ground, body, spring, damping])

        sys = mtkcompile(model)

        prob = ODEProblem(
            sys, [], (0, 20.0); initialization_eqs, fully_determined = true)
        sol = solve(prob; abstol = 1e-9, reltol = 1e-9)

        return sol
    end

    solv = simplify_and_solve(
        dv, sv, bv, gv; initialization_eqs = [bv.s ~ 3, bv.v ~ 1, sv.delta_s ~ 1])
    solp = simplify_and_solve(dp, sp, bp, gp)

    @test solv[bv.v][1] == 1.0
    @test solv[bv.v][end]≈0.0 atol=1e-4

    @test solp[bp.v][1] == 1.0
    @test solp[bp.v][end]≈0.0 atol=1e-4
end

@testset "driven spring damper mass" begin
    @named dv = TV.Damper(d = 1)
    @named dp = TP.Damper(d = 1)

    @named sv = TV.Spring(k = 1)
    @named sp = TP.Spring(k = 1, l = 1)

    @named bv = TV.Mass(m = 1)
    @named bp = TP.Mass(m = 1, v = 1, s = 3)

    @named gv = TV.Fixed()
    @named gp = TP.Fixed(s_0 = 1)

    @named fv = TV.Force()
    @named fp = TP.Force()

    @named source = Sine(frequency = 3, amplitude = 2)

    function TestSystem(damping, spring, body, ground, f, source)
        eqs = [connect(f.f, source.output)
               connect(f.flange, body.flange)
               connect(spring.flange_a, body.flange, damping.flange_a)
               connect(spring.flange_b, damping.flange_b, ground.flange)]

        @named model = System(eqs, t;
            systems = [ground, body, spring, damping, f, source])

        return model
    end

    model = TestSystem(dv, sv, bv, gv, fv, source)
    sys = mtkcompile(model)
    prob = ODEProblem(
        sys, [bv.s => 0, sv.delta_s => 1], (0, 20.0), fully_determined = true)
    solv = solve(prob, Rodas4())

    model = TestSystem(dp, sp, bp, gp, fp, source)
    sys = mtkcompile(model)
    prob = ODEProblem(sys, [], (0, 20.0), fully_determined = true)
    solp = solve(prob, Rodas4())

    for sol in (solv, solp)
        lb, ub = extrema(solv(15:0.05:20, idxs = bv.v).u)
        @test -lb≈ub atol=1e-2
        @test -0.11 < lb < -0.1
    end
end

@testset "sources & sensors" begin
    @testset "Translational" begin
        @testset "PositionSensor & ForceSensor" begin
            function TestSystem(; name)
                systems = @named begin
                    pos = TV.Position()
                    pos_sensor = TV.PositionSensor(; s = 1)
                    force = TV.Force()
                    force_sensor = TV.ForceSensor()

                    spring = TV.Spring(; k = 1000)

                    src1 = Sine(frequency = 100, amplitude = 2)
                    src2 = Sine(frequency = 100, amplitude = -1)

                    pos_value = RealInput()
                    force_output = RealOutput()
                end

                eqs = [connect(pos.s, src1.output)
                       connect(force.f, src2.output)
                       connect(pos.flange, force_sensor.flange_a)
                       connect(force_sensor.flange_b, spring.flange_a)
                       connect(spring.flange_b, force.flange, pos_sensor.flange)
                       connect(pos_value, pos_sensor.output)
                       connect(force_output, force_sensor.output)]

                System(eqs, t, [], []; name, systems)
            end

            @named system = TestSystem()
            s = complete(system)
            sys = mtkcompile(system)
            prob = ODEProblem(sys, [], (0, 1 / 400))
            sol = solve(prob, Rosenbrock23())

            delta_s = 1 / 1000
            s_b = 2 - delta_s + 1

            @test sol[s.pos_value.u][end]≈1.0 atol=1e-3
            @test all(sol[s.spring.flange_a.f] .== sol[s.force_output.u])
        end

        @testset "AccelerationSensor" begin
            @named acc = TV.AccelerationSensor()
            m = 4
            @named mass = TV.Mass(m = m)
            @named force = TV.Force()
            @named source = Sine(frequency = 2, amplitude = 1)
            @named acc_output = RealOutput()
            eqs = [
                connect(force.f, source.output),
                connect(force.flange, mass.flange),
                connect(acc.flange, mass.flange),
                connect(acc_output, acc.output)
            ]
            @named sys = System(
                eqs, t, [], []; systems = [force, source, mass, acc, acc_output])
            s = complete(mtkcompile(sys))
            prob = ODEProblem(s, [mass.s => 0], (0.0, pi))
            sol = solve(prob, Tsit5())
            @test sol[sys.acc_output.u] ≈ (sol[sys.mass.f] ./ m)
        end
    end

    @testset "TranslationalPosition" begin
        @testset "PositionSensor & ForceSensor" begin
            function mass_spring(; name)
                systems = @named begin
                    fixed = TP.Fixed()
                    spring = TP.Spring(; k = 10.0, l = 1.0)
                    mass = TP.Mass(; m = 100.0, s = 2.0, v = 0.0)
                    pos_sensor = TP.PositionSensor()
                    force_sensor = TP.ForceSensor()
                    pos_value = RealOutput()
                    force_value = RealOutput()
                end
                eqs = [
                    connect(fixed.flange, force_sensor.flange_a),
                    connect(force_sensor.flange_b, spring.flange_a),
                    connect(spring.flange_b, mass.flange, pos_sensor.flange),
                    connect(pos_sensor.output, pos_value),
                    connect(force_sensor.output, force_value)
                ]
                System(eqs, t, [], []; name, systems)
            end

            @named model = mass_spring()
            sys = mtkcompile(model)

            prob = ODEProblem(sys, [], (0.0, 1.0), fully_determined = true)
            sol = solve(prob, Tsit5())

            @test all(sol[sys.spring.flange_a.f] .== sol[sys.force_value.u])
            @test all(sol[sys.mass.s] .== sol[sys.pos_value.u])
        end

        @testset "AccelerationSensor" begin
            @named acc = TP.AccelerationSensor()
            m = 4
            @named mass = TP.Mass(m = m)
            @named force = TP.Force()
            @named source = Sine(frequency = 2, amplitude = 1)
            @named acc_output = RealOutput()
            eqs = [
                connect(force.f, source.output),
                connect(force.flange, mass.flange),
                connect(acc.flange, mass.flange),
                connect(acc_output, acc.output)
            ]
            @named sys = System(
                eqs, t, [], []; systems = [force, source, mass, acc, acc_output])
            s = complete(mtkcompile(sys))
            prob = ODEProblem(s, [], (0.0, pi), fully_determined = true)
            sol = solve(prob, Tsit5())
            @test sol[sys.acc_output.u] ≈ (sol[sys.mass.f] ./ m)
        end
    end
end
