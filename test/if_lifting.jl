using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D, IfLifting, no_if_lift

@testset "Simple `abs(x)`" begin
    @mtkmodel SimpleAbs begin
        @variables begin
            x(t)
            y(t)
        end
        @equations begin
            D(x) ~ abs(y)
            y ~ sin(t)
        end
    end
    @named sys = SimpleAbs()
    ss1 = structural_simplify(sys)
    @test length(equations(ss1)) == 1
    ss2 = structural_simplify(sys, additional_passes = [IfLifting])
    @test length(equations(ss2)) == 1
    @test length(parameters(ss2)) == 1
    @test operation(only(equations(ss2)).rhs) === ifelse

    discvar = only(parameters(ss2))
    prob1 = ODEProblem(ss1, [ss1.x => 0.0], (0.0, 5.0))
    sol1 = solve(prob1, Tsit5())
    prob2 = ODEProblem(ss2, [ss2.x => 0.0], (0.0, 5.0))
    sol2 = solve(prob2, Tsit5())
    @test count(isapprox(pi), sol2.t) == 2
    @test any(isapprox(pi), sol2.discretes[1].t)
    @test !sol2[discvar][1]
    @test sol2[discvar][end]

    _t = pi + 1.0
    # x(t) = 1 - cos(t) in [0, pi)
    # x(t) = 3 + cos(t) in [pi, 2pi)
    _trueval = 3 + cos(_t)
    @test !isapprox(sol1(_t)[1], _trueval; rtol = 1e-3)
    @test isapprox(sol2(_t)[1], _trueval; rtol = 1e-3)
end

@testset "Big test case" begin
    @mtkmodel BigModel begin
        @variables begin
            x(t)
            y(t)
            z(t)
            c(t)::Bool
            w(t)
            q(t)
            r(t)
        end
        @parameters begin
            p
        end
        @equations begin
            # ifelse, max, min
            D(x) ~ ifelse(c, max(x, y), min(x, y))
            # discrete observed
            c ~ x <= y
            # observed should also get if-lifting
            y ~ abs(sin(t))
            # should be ignored
            D(z) ~ no_if_lift(ifelse(x < y, x, y))
            # ignore time-independent ifelse
            D(w) ~ ifelse(p < 3, 1.0, 2.0)
            # all the boolean operators
            D(q) ~ ifelse((x < 1) & ((y < 0.5) | ifelse(y > 0.8, c, !c)), 1.0, 2.0)
            # don't touch time-independent condition, but modify time-dependent branches
            D(r) ~ ifelse(p < 2, abs(x), max(y, 0.9))
        end
    end

    @named sys = BigModel()
    ss = structural_simplify(sys, additional_passes = [IfLifting])

    ps = parameters(ss)
    @test length(ps) == 9
    eqs = equations(ss)
    obs = observed(ss)

    @testset "no_if_lift is untouched" begin
        idx = findfirst(eq -> isequal(eq.lhs, D(ss.z)), eqs)
        eq = eqs[idx]
        @test isequal(eq.rhs, no_if_lift(ifelse(ss.x < ss.y, ss.x, ss.y)))
    end
    @testset "time-independent ifelse is untouched" begin
        idx = findfirst(eq -> isequal(eq.lhs, D(ss.w)), eqs)
        eq = eqs[idx]
        @test operation(arguments(eq.rhs)[1]) === Base.:<
    end
    @testset "time-dependent branch of time-independent condition is modified" begin
        idx = findfirst(eq -> isequal(eq.lhs, D(ss.r)), eqs)
        eq = eqs[idx]
        @test operation(eq.rhs) === ifelse
        args = arguments(eq.rhs)
        @test operation(args[1]) == Base.:<
        @test operation(args[2]) === ifelse
        condvars = ModelingToolkit.vars(arguments(args[2])[1])
        @test length(condvars) == 1 && any(isequal(only(condvars)), ps)
        @test operation(args[3]) === ifelse
        condvars = ModelingToolkit.vars(arguments(args[3])[1])
        @test length(condvars) == 1 && any(isequal(only(condvars)), ps)
    end
    @testset "Observed variables are modified" begin
        idx = findfirst(eq -> isequal(eq.lhs, ss.c), obs)
        eq = obs[idx]
        @test operation(eq.rhs) === Base.:! && any(isequal(only(arguments(eq.rhs))), ps)
        idx = findfirst(eq -> isequal(eq.lhs, ss.y), obs)
        eq = obs[idx]
        @test operation(eq.rhs) === ifelse
    end
end

@testset "`@mtkbuild` macro accepts `additional_passes`" begin
    @mtkmodel SimpleAbs begin
        @variables begin
            x(t)
            y(t)
        end
        @equations begin
            D(x) ~ abs(y)
            y ~ sin(t)
        end
    end
    @test_nowarn @mtkbuild sys=SimpleAbs() additional_passes=[IfLifting]
end
