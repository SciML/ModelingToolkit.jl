using ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkit: t_nounits as t, D_nounits as D

using ModelingToolkitStandardLibrary.Blocks: Sine
using ModelingToolkitStandardLibrary.Mechanical.TranslationalModelica: Damper, Spring, Mass,
                                                                       Fixed, Force,
                                                                       SpringDamper,
                                                                       Position

@testset "spring damper mass fixed" begin
    @mtkmodel SpringDamperMassFixed begin
        @components begin
            damper = Damper(; d = 1)
            spring = Spring(; c = 1, s_rel0 = 1)
            mass = Mass(; m = 1, v = 1, s = 0)
            fixed = Fixed(s0 = 1)
        end
        @equations begin
            connect(spring.flange_a, mass.flange_a, damper.flange_a)
            connect(spring.flange_b, damper.flange_b, fixed.flange)
        end
    end

    @mtkcompile sys = SpringDamperMassFixed()

    prob = ODEProblem(sys, [], (0, 20.0))
    sol = solve(prob, ImplicitMidpoint(), dt = 0.01)

    @test sol[sys.mass.v][1] == 1.0
    @test sol[sys.mass.v][end]≈0.0 atol=1e-4
end

@testset "driven spring damper mass" begin
    @mtkmodel DrivenSpringDamperMass begin
        @components begin
            damper = Damper(; d = 1)
            spring = Spring(; c = 1, s_rel0 = 1)
            mass = Mass(; m = 1, v = 1, s = 0)
            fixed = Fixed(; s0 = 1)
            force = Force()
            source = Sine(frequency = 3, amplitude = 2)
        end

        @equations begin
            connect(force.f, source.output)
            connect(force.flange, mass.flange_a)
            connect(spring.flange_a, mass.flange_b, damper.flange_a)
            connect(spring.flange_b, damper.flange_b, fixed.flange)
        end
    end

    @mtkcompile sys = DrivenSpringDamperMass()

    prob = ODEProblem(sys, [], (0, 20.0))
    sol = solve(prob, Rodas4())

    lb, ub = extrema(sol(15:0.05:20, idxs = sys.mass.v).u)
    @test -lb≈ub atol=1e-2
    @test -0.11 < lb < -0.1
end

@testset "driven SpringDamper mass" begin
    @mtkmodel DrivenSpringDamperMass2 begin
        @components begin
            springdamper = SpringDamper(; d = 1, c = 1, s_rel0 = 1)
            mass = Mass(; m = 1, v = 1, s = 0)
            fixed = Fixed(; s0 = 1)
            force = Force()
            source = Sine(frequency = 3, amplitude = 2)
        end

        @equations begin
            connect(force.f, source.output)
            connect(force.flange, mass.flange_a)
            connect(springdamper.flange_a, mass.flange_b)
            connect(springdamper.flange_b, fixed.flange)
        end
    end

    @mtkcompile sys = DrivenSpringDamperMass2()

    prob = ODEProblem(sys, [], (0, 20.0))
    sol = solve(prob, Rodas4())

    lb, ub = extrema(sol(15:0.05:20, idxs = sys.mass.v).u)
    @test -lb≈ub atol=1e-2
    @test -0.11 < lb < -0.1
end

@testset "Position source" begin
    @mtkmodel TestPositionSource begin
        @components begin
            p1 = Position(exact = true)
            source = Sine(frequency = 3, amplitude = 2)
            mass = Mass(m = 1, v = 1, s = 0)
        end

        @equations begin
            connect(source.output, p1.s_ref)
            connect(p1.flange, mass.flange_a)
        end
    end

    @mtkcompile sys = TestPositionSource()
    prob = ODEProblem(sys, [], (0, 2pi))
    sol = solve(prob, Rodas4())
    tv = 0:0.1:(2pi)
    @test sol(tv, idxs = sys.mass.s).u≈@.(2sin(2pi * tv * 3)) atol=1e-2
end
