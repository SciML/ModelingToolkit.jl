using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit: t_nounits as t
using ModelingToolkitStandardLibrary.Blocks: _clamp, _dead_zone
using OrdinaryDiffEq: ReturnCode.Success

@testset "Limiter" begin
    @testset "Constant" begin
        @named c = Constant(; k = 1)
        @named int = Integrator(; k = 1)
        @named sat = Limiter(; y_min = -0.6, y_max = 0.8)
        @named model = System(
            [
                connect(c.output, int.input),
                connect(int.output, sat.input)
            ],
            t,
            systems = [int, c, sat])
        sys = mtkcompile(model)
        prob = ODEProblem(sys, [int.x => 1.0], (0.0, 1.0))

        sol = solve(prob, Rodas4())
        @test SciMLBase.successful_retcode(sol)
        @test sol[int.output.u][end] ≈ 2
        @test sol[sat.output.u][end] ≈ 0.8
    end

    @testset "Sine" begin
        y_min, y_max = -0.3, 0.5
        @named source = Sine(; frequency = 1 / 2)
        @named lim = Limiter(; y_max = y_max, y_min = y_min)
        @named int = Integrator(; k = 1)
        @named iosys = System(
            [
                connect(source.output, lim.input),
                connect(lim.output, int.input)
            ],
            t,
            systems = [source, lim, int])
        sys = mtkcompile(iosys)

        prob = ODEProblem(sys, unknowns(sys) .=> 0.0, (0.0, 10.0))

        sol = solve(prob, Rodas4())
        @test SciMLBase.successful_retcode(sol)
        @test all(abs.(sol[lim.output.u]) .<= 0.5)
        @test all(isapprox.(sol[lim.output.u], _clamp.(sol[source.output.u], y_min, y_max),
            atol = 1e-2))

        # Plots.plot(sol; vars=[source.output.u, lim.output.u])
        # Plots.scatter(sol[source.output.u], sol[lim.output.u])
        # Plots.scatter!(sol[source.output.u], _clamp.(sol[source.output.u], y_min, y_max))
    end
end

@testset "DeadZone" begin
    @testset "Constant" begin
        @named c = Constant(; k = 1)
        @named int = Integrator(; k = 1)
        @named dz = DeadZone(; u_min = -2, u_max = 1)
        @named model = System(
            [
                connect(c.output, int.input),
                connect(int.output, dz.input)
            ],
            t,
            systems = [int, c, dz])
        sys = mtkcompile(model)
        prob = ODEProblem(sys, [int.x => 1.0], (0.0, 1.0))
        sol = solve(prob, Rodas4())

        @test SciMLBase.successful_retcode(sol)
        @test all(sol[int.output.u][end] .≈ 2)
    end

    @testset "Sine" begin
        u_min, u_max = -2, 1
        @named source = Sine(; amplitude = 3, frequency = 1 / 2)
        @named dz = DeadZone(; u_min = u_min, u_max = u_max)
        @named int = Integrator(; k = 1)
        @named model = System(
            [
                connect(source.output, dz.input),
                connect(dz.output, int.input)
            ],
            t,
            systems = [int, source, dz])
        sys = mtkcompile(model)
        prob = ODEProblem(sys, [int.x => 1.0], (0.0, 10.0))
        sol = solve(prob, Rodas4())

        @test SciMLBase.successful_retcode(sol)
        @test all(sol[dz.output.u] .<= 2)
        @test all(sol[dz.output.u] .>= -1)
        @test all(isapprox.(sol[dz.output.u],
            _dead_zone.(sol[source.output.u], u_min, u_max), atol = 1e-2))

        # Plots.plot(sol; vars=[source.output.u, dz.output.u])
        # Plots.scatter(sol[source.output.u], sol[dz.output.u])
        # Plots.scatter!(sol[source.output.u], _dead_zone.(sol[source.output.u], u_min, u_max))
    end
end

@testset "SlewRateLimiter" begin
    @named source = Sine(; frequency = 1 / 2)
    @named rl = SlewRateLimiter(; rising = 1, falling = -1, Td = 0.001, y_start = -1 / 3)
    @named iosys = System([
            connect(source.output, rl.input)
        ],
        t,
        systems = [source, rl])
    sys = mtkcompile(iosys)

    prob = ODEProblem(sys, Pair[], (0.0, 10.0))

    tS = 0.01
    sol = solve(prob, Rodas4(), saveat = tS, abstol = 1e-10, reltol = 1e-10)
    @test SciMLBase.successful_retcode(sol)
    @test all(abs.(sol[rl.output.u]) .<= 0.51)
    @test all(-1 - 1e-5 .<= diff(sol[rl.output.u]) ./ tS .<= 1 + 1e-5) # just an approximation
end
