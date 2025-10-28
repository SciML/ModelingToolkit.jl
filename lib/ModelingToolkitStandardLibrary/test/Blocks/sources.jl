using ModelingToolkit, ModelingToolkitStandardLibrary, OrdinaryDiffEq
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitStandardLibrary.Blocks: smooth_sin, smooth_cos, smooth_damped_sin,
                                             smooth_square, smooth_step, smooth_ramp,
                                             smooth_triangular, triangular, square
using OrdinaryDiffEq: ReturnCode.Success
using DataInterpolations
using DataFrames
using SymbolicIndexingInterface
using SciMLStructures: SciMLStructures, Tunable
using ForwardDiff
using ADTypes

@testset "Constant" begin
    @named src = Constant(k = 2)
    @named int = Integrator()
    @named iosys = System([
            connect(src.output, int.input)
        ],
        t,
        systems = [int, src])
    sys = mtkcompile(iosys)

    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 10.0))

    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test sol[src.output.u][end]≈2 atol=1e-3
end

@testset "TimeVaryingFunction" begin
    f(t) = t^2 + 1
    vars = @variables y(t) dy(t) ddy(t)
    @named src = TimeVaryingFunction(f)
    @named int = Integrator()
    @named iosys = System(
        [y ~ src.output.u
         D(y) ~ dy
         D(dy) ~ ddy
         connect(src.output, int.input)],
        t,
        systems = [int, src])
    sys = mtkcompile(iosys)

    prob = ODEProblem(sys, unknowns(sys) .=> 0.0, (0.0, 10.0))

    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test sol[src.output.u]≈f.(sol.t) atol=1e-3
    @test sol[int.output.u][end]≈1 / 3 * 10^3 + 10 atol=1e-3 # closed-form solution to integral
end

@testset "Sine" begin
    function sine(t, frequency, amplitude, phase, offset, start_time)
        offset + ifelse(t < start_time, 0,
            amplitude * sin(2 * pi * frequency * (t - start_time) + phase))
    end

    frequency = 1
    amplitude = 2
    phase = 0
    offset = 1
    start_time = 2
    δ = 1e-5
    @named int = Integrator()

    @named src = Sine(frequency = frequency, amplitude = amplitude, phase = phase,
        offset = offset, start_time = start_time)
    @named iosys = System([
            connect(src.output, int.input)
        ],
        t,
        systems = [int, src])
    sys = mtkcompile(iosys)

    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 10.0))

    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test sol[src.output.u]≈sine.(sol.t, frequency, amplitude, phase, offset, start_time) atol=1e-3

    @named smooth_src = Sine(frequency = frequency,
        amplitude = amplitude,
        phase = phase,
        offset = offset,
        start_time = start_time,
        smooth = true)
    @named smooth_iosys = System([
            connect(smooth_src.output, int.input)
        ],
        t,
        systems = [int, smooth_src])

    smooth_sys = mtkcompile(smooth_iosys)
    smooth_prob = ODEProblem(smooth_sys, Pair[int.x => 0.0], (0.0, 10.0))
    smooth_sol = solve(smooth_prob, Rodas4())

    @test sol.retcode == Success
    @test smooth_sol[smooth_src.output.u]≈smooth_sin.(
        smooth_sol.t, δ, frequency, amplitude,
        phase, offset, start_time) atol=1e-3
end

@testset "Cosine" begin
    function cosine(t, frequency, amplitude, phase, offset, start_time)
        offset + ifelse(t < start_time, 0,
            amplitude * cos(2 * pi * frequency * (t - start_time) + phase))
    end

    frequency = 1
    amplitude = 2
    phase = 0
    offset = 1
    start_time = 2
    δ = 1e-5
    @named int = Integrator()

    @named src = Cosine(frequency = frequency,
        amplitude = amplitude,
        phase = phase,
        offset = offset,
        start_time = start_time,
        smooth = false)
    @named iosys = System([
            connect(src.output, int.input)
        ],
        t,
        systems = [int, src])

    sys = mtkcompile(iosys)
    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 10.0))
    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test sol[src.output.u]≈cosine.(sol.t, frequency, amplitude, phase, offset, start_time) atol=1e-3

    @named smooth_src = Cosine(frequency = frequency,
        amplitude = amplitude,
        phase = phase,
        offset = offset,
        start_time = start_time,
        smooth = true)
    @named smooth_iosys = System([
            connect(smooth_src.output, int.input)
        ],
        t,
        systems = [int, smooth_src])

    smooth_sys = mtkcompile(smooth_iosys)
    smooth_prob = ODEProblem(smooth_sys, Pair[int.x => 0.0], (0.0, 10.0))
    smooth_sol = solve(smooth_prob, Rodas4())

    @test smooth_sol.retcode == Success
    @test smooth_sol[smooth_src.output.u]≈smooth_cos.(
        smooth_sol.t, δ, frequency, amplitude,
        phase, offset, start_time) atol=1e-3
end

@testset "ContinuousClock" begin
    cont_clock(t, offset, start_time) = offset + ifelse(t < start_time, 0, t - start_time)

    offset, start_time = 1, 0

    @named src = ContinuousClock(offset = offset, start_time = start_time)
    @named int = Integrator()
    @named iosys = System([
            connect(src.output, int.input)
        ],
        t,
        systems = [int, src])
    sys = mtkcompile(iosys)

    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 10.0))

    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test sol[src.output.u]≈cont_clock.(sol.t, offset, start_time) atol=1e-3
end

@testset "Ramp" begin
    function ramp(t, offset, height, duration, start_time)
        offset + ifelse(t < start_time, 0,
            ifelse(t < (start_time + duration), (t - start_time) * height / duration,
                height))
    end

    offset, height, duration, start_time, δ = 1, 2, 2, 0, 1e-5
    @named int = Integrator()

    @named src = Ramp(offset = offset, height = height, duration = duration,
        start_time = start_time)
    @named iosys = System([
            connect(src.output, int.input)
        ],
        t,
        systems = [int, src])

    sys = mtkcompile(iosys)
    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 10.0))
    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test sol[src.output.u]≈ramp.(sol.t, offset, height, duration, start_time) atol=1e-3

    start_time = 2
    @named smooth_src = Ramp(offset = offset, height = height, duration = duration,
        start_time = start_time, smooth = true)
    @named smooth_iosys = System([
            connect(smooth_src.output, int.input)
        ],
        t,
        systems = [int, smooth_src])

    smooth_sys = mtkcompile(smooth_iosys)
    smooth_prob = ODEProblem(smooth_sys, Pair[int.x => 0.0], (0.0, 10.0))
    smooth_sol = solve(smooth_prob, Rodas4())

    @test smooth_sol.retcode == Success
    @test smooth_sol[smooth_src.output.u]≈smooth_ramp.(smooth_sol.t, δ, height, duration,
        offset, start_time) atol=1e-3
end

@testset "Step" begin
    step(t, offset, height, start_time) = offset + ifelse(t < start_time, 0, height)

    offset, height, start_time, δ = 1, 2, 5, 1e-5
    @named int = Integrator()

    @named src = Step(offset = offset, height = height, start_time = start_time,
        smooth = false)
    @named iosys = System([
            connect(src.output, int.input)
        ],
        t,
        systems = [int, src])

    sys = mtkcompile(iosys)
    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 10.0))
    sol = solve(prob, Rodas4())

    @test sol.retcode == Success
    @test sol[src.output.u]≈step.(sol.t, offset, height, start_time) atol=1e-2
    @test sol(start_time, idxs = src.output.u) == height + offset # Test that the step is applied at the start time

    # test with duration
    duration = 1.2
    @named src = Step(offset = offset, height = height, start_time = start_time,
        duration = duration, smooth = false)
    @named iosys = System([
            connect(src.output, int.input)
        ],
        t,
        systems = [int, src])

    sys = mtkcompile(iosys)
    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 10.0))
    sol = solve(prob, Rodas4(), dtmax = 0.1) # set dtmax to prevent the solver from overstepping the entire step disturbance

    @test sol.retcode == Success
    @test sol[src.output.u]≈step.(sol.t, offset, height, start_time) -
                            step.(sol.t, 0, height, start_time + duration) atol=1e-2

    @named smooth_src = Step(offset = offset, height = height, start_time = start_time,
        smooth = true)
    @named smooth_iosys = System([
            connect(smooth_src.output, int.input)
        ],
        t,
        systems = [int, smooth_src])

    smooth_sys = mtkcompile(smooth_iosys)
    smooth_prob = ODEProblem(smooth_sys, Pair[int.x => 0.0], (0.0, 10.0))
    smooth_sol = solve(smooth_prob, Rodas4(), dtmax = 0.1) # set dtmax to prevent the solver from overstepping the entire step disturbance)

    @test smooth_sol.retcode == Success
    @test smooth_sol[smooth_src.output.u] ≈
          smooth_step.(smooth_sol.t, δ, height, offset, start_time)

    # with duration
    @named smooth_src = Step(offset = offset, height = height, start_time = start_time,
        smooth = true, duration = duration)
    @named smooth_iosys = System([
            connect(smooth_src.output, int.input)
        ],
        t,
        systems = [int, smooth_src])

    smooth_sys = mtkcompile(smooth_iosys)
    smooth_prob = ODEProblem(smooth_sys, Pair[int.x => 0.0], (0.0, 10.0))
    smooth_sol = solve(smooth_prob, Rodas4())

    @test smooth_sol.retcode == Success
    @test smooth_sol[smooth_src.output.u] ≈
          smooth_step.(smooth_sol.t, δ, height, offset, start_time) -
          smooth_step.(smooth_sol.t, δ, height, 0, start_time + duration)
end

@testset "Square" begin
    frequency = 1
    amplitude = 2
    offset = 1
    start_time = 2.5
    δ = 1e-5
    @named int = Integrator()

    @named src = Square(frequency = frequency, amplitude = amplitude,
        offset = offset, start_time = start_time)
    @named iosys = System([
            connect(src.output, int.input)
        ],
        t,
        systems = [int, src])

    sys = mtkcompile(iosys)
    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 10.0))
    sol = solve(prob, Rodas4())

    @test sol.retcode == Success
    @test sol[src.output.u]≈square.(sol.t, frequency, amplitude, offset, start_time) atol=1e-3

    @named smooth_src = Square(frequency = frequency, amplitude = amplitude,
        offset = offset, start_time = start_time, smooth = true)
    @named smooth_iosys = System([
            connect(smooth_src.output, int.input)
        ],
        t,
        systems = [int, smooth_src])

    smooth_sys = mtkcompile(smooth_iosys)
    smooth_prob = ODEProblem(smooth_sys, Pair[int.x => 0.0], (0.0, 10.0))
    smooth_sol = solve(smooth_prob, Rodas4())

    @test smooth_sol.retcode == Success
    @test smooth_sol[smooth_src.output.u]≈smooth_square.(smooth_sol.t, δ, frequency,
        amplitude, offset, start_time) atol=1e-3
end

@testset "Triangular" begin
    frequency = 5
    amplitude = 1
    offset = 2
    start_time = 1
    δ = 1e-5
    @named int = Integrator()

    @named src = Triangular(frequency = frequency, amplitude = amplitude,
        offset = offset, start_time = start_time)
    @named iosys = System([
            connect(src.output, int.input)
        ],
        t,
        systems = [int, src])

    sys = mtkcompile(iosys)
    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 4.0))
    sol = solve(prob, Rodas4(), saveat = 0.01)

    @test sol.retcode == Success
    @test sol[src.output.u]≈triangular.(sol.t, frequency, amplitude, offset, start_time) atol=1e-3

    @named smooth_src = Triangular(frequency = frequency, amplitude = amplitude,
        offset = offset, start_time = start_time, smooth = true)
    @named smooth_iosys = System([
            connect(smooth_src.output, int.input)
        ],
        t,
        systems = [int, smooth_src])

    smooth_sys = mtkcompile(smooth_iosys)
    smooth_prob = ODEProblem(smooth_sys, Pair[int.x => 0.0], (0.0, 4.0))
    smooth_sol = solve(smooth_prob, Rodas4(), saveat = 0.01)

    @test smooth_sol.retcode == Success
    @test smooth_sol[smooth_src.output.u]≈smooth_triangular.(smooth_sol.t, δ, frequency,
        amplitude, offset, start_time) atol=1e-3
end

@testset "ExpSine" begin
    function exp_sine(t, amplitude, frequency, damping, phase, start_time)
        offset + ifelse(t < start_time, 0,
            amplitude * exp(-damping * (t - start_time)) *
            sin(2 * pi * frequency * (t - start_time) + phase))
    end

    frequency, amplitude, damping = 3, 2, 0.10
    phase, offset, start_time, δ = 0, 0, 0, 1e-5
    @named src = ExpSine(frequency = frequency, amplitude = amplitude, damping = damping,
        phase = phase, offset = offset, start_time = start_time)
    @named int = Integrator()
    @named iosys = System([
            connect(src.output, int.input)
        ],
        t,
        systems = [int, src])
    sys = mtkcompile(iosys)
    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 10.0))
    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test sol[src.output.u]≈exp_sine.(sol.t, amplitude, frequency, damping, phase,
        start_time) atol=1e-3

    offset, start_time = 1, 2
    @named smooth_src = ExpSine(frequency = frequency, amplitude = amplitude,
        damping = damping, phase = phase, offset = offset,
        start_time = start_time, smooth = true)
    @named smooth_iosys = System([
            connect(smooth_src.output, int.input)
        ],
        t,
        systems = [int, smooth_src])
    smooth_sys = mtkcompile(smooth_iosys)
    smooth_prob = ODEProblem(smooth_sys, Pair[int.x => 0.0], (0.0, 10.0))
    smooth_sol = solve(smooth_prob, Rodas4())

    @test smooth_sol.retcode == Success
    @test smooth_sol[smooth_src.output.u]≈smooth_damped_sin.(smooth_sol.t, δ, frequency,
        amplitude, damping, phase,
        offset, start_time) atol=1e-3
end

@testset "SampledData" begin
    dt = 4e-4
    t_end = 10.0
    time = 0:dt:t_end
    x = @. time^2 + 1.0

    @testset "using Parameter type" begin
        vars = @variables y(t) dy(t) ddy(t)
        @named src = SampledData(Float64)
        @named int = Integrator()
        @named iosys = System(
            [y ~ src.output.u
             D(y) ~ dy
             D(dy) ~ ddy
             connect(src.output, int.input)],
            t,
            systems = [int, src])
        sys = mtkcompile(iosys)
        s = complete(iosys)
        prob = ODEProblem(sys,
            [s.src.buffer => Parameter(x, dt)],
            (0.0, t_end);
            tofloat = false)
        # prob = remake(prob; p = Parameter.(prob.p)) #<-- no longer needed with ModelingToolkit.jl PR #2231

        sol = solve(prob, Rodas4())
        @test sol.retcode == Success
        @test sol[src.output.u][1] == 1.0 #check correct initial condition

        @test sol(time)[src.output.u]≈x atol=1e-3
        @test sol[int.output.u][end]≈1 / 3 * 10^3 + 10.0 atol=1e-3 # closed-form solution to integral
        @test sol[dy][end]≈2 * time[end] atol=1e-3
        @test sol[ddy][end]≈2 atol=1e-3
    end

    @testset "using Vector Based" begin
        vars = @variables y(t) dy(t) ddy(t)
        @named src = SampledData(dt)
        @named int = Integrator()
        @named iosys = System(
            [y ~ src.output.u
             D(y) ~ dy
             D(dy) ~ ddy
             connect(src.output, int.input)],
            t,
            systems = [int, src])
        sys = mtkcompile(iosys)
        s = complete(iosys)
        prob = ODEProblem(sys,
            [s.src.buffer => x, s.src.sample_time => dt],
            (0.0, t_end);
            tofloat = false)

        sol = solve(prob, Rodas4())
        @test sol.retcode == Success
        @test sol[src.output.u][1] == 1.0 #check correct initial condition

        @test sol(time)[src.output.u]≈x atol=1e-3
        @test sol[int.output.u][end]≈1 / 3 * 10^3 + 10.0 atol=1e-3 # closed-form solution to integral
        @test sol[dy][end]≈2 * time[end] atol=1e-3
        @test sol[ddy][end]≈2 atol=1e-3
    end
end

@testset "Interpolation" begin
    @variables y(t) = 0
    u = rand(15)
    x = 0:14.0

    @named i = Interpolation(LinearInterpolation, u, x)
    eqs = [i.input.u ~ t, D(y) ~ i.output.u]

    @named model = System(eqs, t, systems = [i])
    sys = mtkcompile(model)

    prob = ODEProblem{true, SciMLBase.FullSpecialize}(sys, [], (0.0, 4))
    sol = solve(prob, Tsit5())

    @test SciMLBase.successful_retcode(sol)
end

@testset "Interpolation in model macro" begin
    function MassSpringDamper(; name)
        @named input = RealInput()
        @variables f(t) x(t)=0 dx(t)=0 ddx(t)
        @parameters m=10 k=1000 d=1

        eqs = [f ~ input.u
               ddx * 10 ~ k * x + d * dx + f
               D(x) ~ dx
               D(dx) ~ ddx]

        System(eqs, t; name, systems = [input])
    end

    table_data = [1.0, 2.0, 3.0]
    table_bkp = [0.0, 0.5, 1.0]
    itp = LinearInterpolation(table_data, table_bkp)

    @mtkmodel model_with_lut begin
        @components begin
            src = Interpolation(itp)
            clk = ContinuousClock()
            model = MassSpringDamper()
        end
        @equations begin
            connect(src.input, clk.output)
            connect(src.output, model.input)
        end
    end;
    @mtkcompile sys = model_with_lut()

    prob = ODEProblem(sys, [], (0.0, 1))
    sol = solve(prob, Tsit5())

    @test SciMLBase.successful_retcode(sol)
end

@testset "ParametrizedInterpolation" begin
    @variables y(t) = 0
    u = rand(15)
    x = 0:14.0

    @testset "LinearInterpolation" begin
        @named i = ParametrizedInterpolation(LinearInterpolation, u, x)
        eqs = [i.input.u ~ t, D(y) ~ i.output.u]

        @named model = System(eqs, t, systems = [i])
        sys = mtkcompile(model)

        prob = ODEProblem{true, SciMLBase.FullSpecialize}(sys, [], (0.0, 4))
        sol = solve(prob, Tsit5())

        @test SciMLBase.successful_retcode(sol)

        prob2 = remake(prob, p = [i.data => ones(15)])
        sol2 = solve(prob2)

        @test SciMLBase.successful_retcode(sol2)
        @test all(only.(sol2.u) .≈ sol2.t) # the solution for y' = 1 is y(t) = t

        set_data! = setp(prob2, i.data)
        set_data!(prob2, zeros(15))
        sol3 = solve(prob2)
        @test SciMLBase.successful_retcode(sol3)
        @test iszero(sol3)

        function loss(x, p)
            prob0, set_data! = p
            ps = parameter_values(prob0)
            arr, repack, alias = SciMLStructures.canonicalize(Tunable(), ps)
            T = promote_type(eltype(x), eltype(arr))
            promoted_ps = SciMLStructures.replace(Tunable(), ps, T.(arr))
            prob = remake(prob0; p = promoted_ps)

            set_data!(prob, x)
            sol = solve(prob)
            sum(abs2.(only.(sol.u) .- sol.t))
        end

        set_data! = setp(prob, i.data)
        of = OptimizationFunction(loss, AutoForwardDiff())
        op = OptimizationProblem(
            of, u, (prob, set_data!), lb = zeros(15), ub = fill(2.0, 15))

        # check that type changing works
        @test length(ForwardDiff.gradient(x -> of(x, (prob, set_data!)), u)) == 15

        @test_skip begin
            r = solve(op, Optimization.LBFGS(), maxiters = 1000)
            @test of(r.u, (prob, set_data!)) < of(u, (prob, set_data!))
        end
    end

    @testset "BSplineInterpolation" begin
        @named i = ParametrizedInterpolation(
            BSplineInterpolation, u, x, 3, :Uniform, :Uniform)
        eqs = [i.input.u ~ t, D(y) ~ i.output.u]

        @named model = System(eqs, t, systems = [i])
        sys = mtkcompile(model)

        prob = ODEProblem(sys, [], (0.0, 4))
        sol = solve(prob)

        @test SciMLBase.successful_retcode(sol)
    end

    @testset "Initialization" begin
        function MassSpringDamper(; name)
            @named input = RealInput()
            vars = @variables f(t) x(t)=0 dx(t) [guess = 0] ddx(t)
            pars = @parameters m=10 k=1000 d=1

            eqs = [f ~ input.u
                   ddx * 10 ~ k * x + d * dx + f
                   D(x) ~ dx
                   D(dx) ~ ddx]

            System(eqs, t, vars, pars; name, systems = [input])
        end

        function MassSpringDamperSystem(data, time; name)
            @named src = ParametrizedInterpolation(LinearInterpolation, data, time)
            @named clk = ContinuousClock()
            @named model = MassSpringDamper()

            eqs = [connect(model.input, src.output)
                   connect(src.input, clk.output)]

            System(eqs, t; name, systems = [src, clk, model])
        end

        function generate_data()
            dt = 4e-4
            time = 0:dt:0.1
            data = sin.(2 * pi * time * 100)

            return DataFrame(; time, data)
        end

        df = generate_data() # example data

        @named system = MassSpringDamperSystem(df.data, df.time)
        sys = mtkcompile(system)
        prob = ODEProblem(sys, [], (0, df.time[end]))
        sol = solve(prob)

        @test SciMLBase.successful_retcode(sol)

        prob2 = remake(prob, p = [sys.src.data => ones(length(df.data))])
        sol2 = solve(prob2)

        @test SciMLBase.successful_retcode(sol2)
    end
end
