using ModelingToolkit, ModelingToolkitStandardLibrary, OrdinaryDiffEq
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit: t_nounits as t
using OrdinaryDiffEq: ReturnCode.Success
using Test

#=
Testing strategy:
The general strategy is to test systems using simple inputs where the solution
is known on closed form. For algebraic systems (without differential variables),
an integrator with a constant input is often used together with the system under test.
=#

@testset "Constant" begin
    @named c = Constant(; k = 1)
    @named int = Integrator(x = 1)
    @named iosys = System(connect(c.output, int.input), t, systems = [int, c])
    sys = mtkcompile(iosys)
    prob = ODEProblem(sys, Pair[], (0.0, 1.0))
    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test all(sol[c.output.u] .≈ 1)
    @test sol[int.output.u][end] .≈ 2 # expected solution
end

@testset "Derivative" begin
    @named source = Sine(; frequency = 1)
    @named int = Integrator(; k = 1)
    @named der = Derivative(; k = 1, T = 0.001)
    @named iosys = System(
        [
            connect(source.output, der.input),
            connect(der.output, int.input)
        ],
        t,
        systems = [int, source, der])
    sys = mtkcompile(iosys)
    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 10.0))
    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test all(isapprox.(sol[source.output.u], sol[int.output.u], atol = 1e-1))
end

@testset "PT1" begin
    pt1_func(t, k, T) = k * (1 - exp(-t / T)) # Known solution to first-order system

    k, T = 1.2, 0.1
    @named c = Constant(; k = 1)
    @named pt1 = FirstOrder(; k = k, T = T)
    @named iosys = System(connect(c.output, pt1.input), t, systems = [pt1, c])
    sys = mtkcompile(iosys)
    prob = ODEProblem(sys, Pair[], (0.0, 100.0))
    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test sol[pt1.output.u]≈pt1_func.(sol.t, k, T) atol=1e-3

    # Test highpass feature
    @named pt1 = FirstOrder(; k = k, T = T, lowpass = false)
    @named iosys = System(connect(c.output, pt1.input), t, systems = [pt1, c])
    sys = mtkcompile(iosys)
    prob = ODEProblem(sys, Pair[], (0.0, 100.0))
    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test sol[pt1.output.u]≈k .- pt1_func.(sol.t, k, T) atol=1e-3
end

@testset "PT2" begin
    # Known solution to second-order system
    function pt2_func(t, k, w, d)
        y = if d == 0
            -k * (-1 + cos(t * w))
        else
            d = complex(d)
            real(k * (1 +
                  (-cosh(sqrt(-1 + d^2) * t * w) -
                   (d * sinh(sqrt(-1 + d^2) * t * w)) / sqrt(-1 + d^2)) / exp(d * t * w)))
        end
    end

    k, w, d = 1.0, 1.0, 0.5
    @named c = Constant(; k = 1)
    @named pt2 = SecondOrder(; k = k, w = w, d = d)
    @named iosys = System(connect(c.output, pt2.input), t, systems = [pt2, c])
    sys = mtkcompile(iosys)
    prob = ODEProblem(sys, [unknowns(sys) .=> 0.0...; pt2.xd => 0.0], (0.0, 100.0))
    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test sol[pt2.output.u]≈pt2_func.(sol.t, k, w, d) atol=1e-3
end

@testset "StateSpace" begin
    A = [0 1; -1 -0.5]
    B = [0, 1]
    C = [0.9 1;]
    D = [0;;]
    @named ss = StateSpace(; A, B, C, D, x = zeros(2))
    @named c = Constant(; k = 1)
    @named model = System([
            connect(c.output, ss.input)
        ],
        t,
        systems = [ss, c])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, Pair[], (0.0, 100.0))
    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    # initial condition
    @test sol[ss.x[1]][1]≈0 atol=1e-3
    @test sol[ss.x[2]][1]≈0 atol=1e-3
    # equilibrium point is at [1, 0]
    @test sol[ss.x[1]][end]≈1 atol=1e-3
    @test sol[ss.x[2]][end]≈0 atol=1e-3

    # non-zero operating point
    u0 = [1] # This causes no effective input to the system since c.k = 1
    y0 = [2]
    @named ss = StateSpace(; A, B, C, D, x = zeros(2), u0, y0)
    @named model = System([
            connect(c.output, ss.input)
        ],
        t,
        systems = [ss, c])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, Pair[], (0.0, 100.0))
    sol = solve(prob, Rodas4())
    @test sol.retcode == Success

    @test sol[ss.x[1]][end ÷ 2]≈0 atol=1e-3 # Test that x did not move
    @test sol[ss.x[1]][end]≈0 atol=1e-3 # Test that x did not move
    @test sol[ss.x[2]][end]≈0 atol=1e-3
    @test sol[ss.output.u[1]][end]≈y0[] atol=1e-3 # Test that the output equals the operating point
end

"""
Second order demo plant
"""
@component function Plant(; name, x = zeros(2))
    @named input = RealInput()
    @named output = RealOutput()
    D = Differential(t)
    sts = @variables x1(t)=x[1] x2(t)=x[2]
    eqs = [D(x1) ~ x2
           D(x2) ~ -x1 - 0.5 * x2 + input.u
           output.u ~ 0.9 * x1 + x2]
    compose(System(eqs, t, sts, []; name), [input, output])
end

@testset "PI" begin
    re_val = 2
    @named ref = Constant(; k = re_val)
    @named pi_controller = PI(k = 1, T = 1)
    @named plant = Plant()
    @named fb = Feedback()
    @named model = System(
        [
            connect(ref.output, fb.input1),
            connect(plant.output, fb.input2),
            connect(fb.output, pi_controller.err_input),
            connect(pi_controller.ctr_output, plant.input)
        ],
        t,
        systems = [pi_controller, plant, ref, fb])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, Pair[], (0.0, 100.0))
    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test all(isapprox.(sol[ref.output.u], re_val, atol = 1e-3))  # check reference
    @test sol[plant.output.u][end]≈re_val atol=1e-3 # zero control error after 100s
end

@testset "PID" begin
    re_val = 2
    @named ref = Constant(; k = re_val)
    @named pid_controller = PID(k = 3, Ti = 0.5, Td = 1 / 100)
    @named plant = Plant()
    @named fb = Feedback()
    @named model = System(
        [
            connect(ref.output, fb.input1),
            connect(plant.output, fb.input2),
            connect(fb.output, pid_controller.err_input),
            connect(pid_controller.ctr_output, plant.input)
        ],
        t,
        systems = [pid_controller, plant, ref, fb])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, Pair[], (0.0, 100.0))
    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test all(isapprox.(sol[ref.output.u], re_val, atol = 1e-3))  # check reference
    @test sol[plant.output.u][end]≈re_val atol=1e-3 # zero control error after 100s

    @testset "PI" begin
        @named pid_controller = PID(k = 3, Ti = 0.5, Td = false)
        @named model = System(
            [
                connect(ref.output, fb.input1),
                connect(plant.output, fb.input2),
                connect(fb.output, pid_controller.err_input),
                connect(pid_controller.ctr_output, plant.input)
            ],
            t,
            systems = [pid_controller, plant, ref, fb])
        sys = mtkcompile(model)
        prob = ODEProblem(sys, Pair[], (0.0, 100.0))
        sol = solve(prob, Rodas4())
        @test sol.retcode == Success
        @test all(isapprox.(sol[ref.output.u], re_val, atol = 1e-3))  # check reference
        @test sol[plant.output.u][end]≈re_val atol=1e-3 # zero control error after 100s
    end

    @testset "PD" begin
        @named pid_controller = PID(k = 10, Ti = false, Td = 1)
        @named model = System(
            [
                connect(ref.output, fb.input1),
                connect(plant.output, fb.input2),
                connect(fb.output, pid_controller.err_input),
                connect(pid_controller.ctr_output, plant.input)
            ],
            t,
            systems = [pid_controller, plant, ref, fb])
        sys = mtkcompile(model)
        prob = ODEProblem(sys, Pair[], (0.0, 100.0))
        sol = solve(prob, Rodas4())
        @test sol.retcode == Success
        @test all(isapprox.(sol[ref.output.u], re_val, atol = 1e-3))  # check reference
        @test sol[plant.output.u][end] > 1 # without I there will be a steady-state error
    end
end

@testset "LimPI" begin
    re_val = 1
    @named ref = Constant(; k = re_val)
    @named pi_controller_lim = LimPI(k = 3,
        T = 0.5,
        u_max = 1.5,
        u_min = -1.5,
        Ta = 0.1)
    @named pi_controller = PI(gainPI.k = 3, T = 0.5)
    @named sat = Limiter(y_max = 1.5, y_min = -1.5)
    @named plant = Plant()
    @named fb = Feedback()

    # without anti-windup measure
    sol = let
        @named model = System(
            [
                connect(ref.output, fb.input1),
                connect(plant.output, fb.input2),
                connect(fb.output, pi_controller.err_input),
                connect(pi_controller.ctr_output, sat.input),
                connect(sat.output, plant.input)
            ],
            t,
            systems = [pi_controller, plant, ref, fb, sat])
        sys = mtkcompile(model)
        prob = ODEProblem(sys, Pair[], (0.0, 20.0))
        sol = solve(prob, Rodas4())
    end

    # with anti-windup measure
    sol_lim = let
        @named model = System(
            [
                connect(ref.output, fb.input1),
                connect(plant.output, fb.input2),
                connect(fb.output, pi_controller_lim.err_input),
                connect(pi_controller_lim.ctr_output, sat.input),
                connect(sat.output, plant.input)
            ],
            t,
            systems = [pi_controller_lim, plant, ref, fb, sat])
        sys = mtkcompile(model)
        prob = ODEProblem(sys, Pair[], (0.0, 20.0))
        sol = solve(prob, Rodas4())
    end

    @test sol.retcode == Success
    @test sol_lim.retcode == ReturnCode.Success
    @test all(isapprox.(sol[ref.output.u], re_val, atol = 1e-3))  # check reference
    @test all(isapprox.(sol_lim[ref.output.u], re_val, atol = 1e-3))  # check reference
    @test sol[plant.output.u][end]≈re_val atol=1e-3 # zero control error after 100s
    @test sol_lim[plant.output.u][end]≈re_val atol=1e-3 # zero control error after 100s
    @test all(-1.5 .<= sol_lim[pi_controller_lim.ctr_output.u] .<= 1.5) # test limit

    # Plots.plot(sol; vars=[plant.output.u]) # without anti-windup measure
    # Plots.plot!(sol_lim; vars=[plant.output.u]) # with anti-windup measure
end

@testset "LimPID" begin
    re_val = 1
    @named ref = Constant(; k = re_val)
    @named pid_controller = LimPID(
        k = 3, Ti = 0.5, Td = 1 / 100, u_max = 1.5, u_min = -1.5,
        Ni = 0.1 / 0.5)
    @named plant = Plant()
    @named model = System(
        [
            connect(ref.output, pid_controller.reference),
            connect(plant.output, pid_controller.measurement),
            connect(pid_controller.ctr_output, plant.input)
        ],
        t,
        systems = [pid_controller, plant, ref])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, Pair[], (0.0, 100.0))
    sol = solve(prob, Rodas4())

    # Plots.plot(sol, vars=[plant.output.u, plant.input.u])
    @test sol.retcode == Success
    @test all(isapprox.(sol[ref.output.u], re_val, atol = 1e-3))  # check reference
    @test sol[plant.output.u][end]≈re_val atol=1e-3 # zero control error after 100s
    @test all(-1.5 .<= sol[pid_controller.ctr_output.u] .<= 1.5) # test limit

    @testset "PI" begin
        @named pid_controller = LimPID(k = 3, Ti = 0.5, Td = false, u_max = 1.5,
            u_min = -1.5, Ni = 0.1 / 0.5)
        @named model = System(
            [
                connect(ref.output, pid_controller.reference),
                connect(plant.output, pid_controller.measurement),
                connect(pid_controller.ctr_output, plant.input)
            ],
            t,
            systems = [pid_controller, plant, ref])
        sys = mtkcompile(model)
        prob = ODEProblem(sys, Pair[], (0.0, 100.0))
        sol = solve(prob, Rodas4())

        # Plots.plot(sol, vars=[plant.output.u, plant.input.u])
        @test sol.retcode == Success
        @test all(isapprox.(sol[ref.output.u], re_val, atol = 1e-3))  # check reference
        @test sol[plant.output.u][end]≈re_val atol=1e-3 # zero control error after 100s
        @test all(-1.5 .<= sol[pid_controller.ctr_output.u] .<= 1.5) # test limit
    end
    @testset "PD" begin
        @named pid_controller = LimPID(k = 10, Ti = false, Td = 1, u_max = 1.5,
            u_min = -1.5)
        @named model = System(
            [
                connect(ref.output, pid_controller.reference),
                connect(plant.output, pid_controller.measurement),
                connect(pid_controller.ctr_output, plant.input)
            ],
            t,
            systems = [pid_controller, plant, ref])
        sys = mtkcompile(model)
        prob = ODEProblem(sys, Pair[], (0.0, 100.0))
        sol = solve(prob, Rodas4())

        # Plots.plot(sol, vars=[plant.output.u, plant.input.u])
        @test sol.retcode == Success
        @test all(isapprox.(sol[ref.output.u], re_val, atol = 1e-3))  # check reference
        @test sol[plant.output.u][end] > 0.5 # without I there will be a steady-state error
        @test all(-1.5 .<= sol[pid_controller.ctr_output.u] .<= 1.5) # test limit
    end
    @testset "set-point weights" begin
        @testset "wp" begin
            @named pid_controller = LimPID(k = 3, Ti = 0.5, Td = 1 / 100, u_max = 1.5,
                u_min = -1.5, Ni = 0.1 / 0.5, wp = 0, wd = 1)
            @named model = System(
                [
                    connect(ref.output, pid_controller.reference),
                    connect(plant.output, pid_controller.measurement),
                    connect(pid_controller.ctr_output, plant.input)
                ],
                t,
                systems = [pid_controller, plant, ref])
            sys = mtkcompile(model)
            prob = ODEProblem(sys, Pair[], (0.0, 100.0))
            sol = solve(prob, Rodas4())

            # Plots.plot(sol, vars=[plant.output.u, plant.input.u])
            @test sol.retcode == Success
            @test all(isapprox.(sol[ref.output.u], re_val, atol = 1e-3))  # check reference
            sol[pid_controller.addP.output.u] == -sol[pid_controller.measurement.u]
            @test sol[plant.output.u][end]≈re_val atol=1e-3 # zero control error after 100s
            @test all(-1.5 .<= sol[pid_controller.ctr_output.u] .<= 1.5) # test limit
        end
        @testset "wd" begin
            @named pid_controller = LimPID(k = 3, Ti = 0.5, Td = 1 / 100, u_max = 1.5,
                u_min = -1.5, Ni = 0.1 / 0.5, wp = 1, wd = 0)
            @named model = System(
                [
                    connect(ref.output, pid_controller.reference),
                    connect(plant.output, pid_controller.measurement),
                    connect(pid_controller.ctr_output, plant.input)
                ],
                t,
                systems = [pid_controller, plant, ref])
            sys = mtkcompile(model)
            prob = ODEProblem(sys, Pair[], (0.0, 100.0))
            sol = solve(prob, Rodas4())

            # Plots.plot(sol, vars=[plant.output.u, plant.input.u])
            @test sol.retcode == Success
            @test all(isapprox.(sol[ref.output.u], re_val, atol = 1e-3))  # check reference
            @test sol[plant.output.u][end]≈re_val atol=1e-3 # zero control error after 100s
            sol[pid_controller.addD.output.u] == -sol[pid_controller.measurement.u]
            @test all(-1.5 .<= sol[pid_controller.ctr_output.u] .<= 1.5) # test limit
        end
    end
    @testset "PI without AWM" begin
        @named pid_controller = LimPID(k = 3, Ti = 0.5, Td = false, u_max = 1.5,
            u_min = -1.5, Ni = Inf)
        @named model = System(
            [
                connect(ref.output, pid_controller.reference),
                connect(plant.output, pid_controller.measurement),
                connect(pid_controller.ctr_output, plant.input)
            ],
            t,
            systems = [pid_controller, plant, ref])
        sys = mtkcompile(model)
        prob = ODEProblem(sys, Pair[], (0.0, 100.0))
        sol = solve(prob, Rodas4())

        # Plots.plot(sol, vars=[plant.output.u, plant.input.u])
        @test sol.retcode == Success
        @test all(isapprox.(sol[ref.output.u], re_val, atol = 1e-3))  # check reference
        @test sol[plant.output.u][end]≈re_val atol=1e-3 # zero control error after 100s
        @test all(-1.5 .<= sol[pid_controller.ctr_output.u] .<= 1.5) # test limit
    end

    @testset "TransferFunction" begin
        pt1_func(t, k, T) = k * (1 - exp(-t / T)) # Known solution to first-order system

        @named c = Constant(; k = 1)
        @named pt1 = TransferFunction(b = [1.2], a = [3.14, 1])
        @named iosys = System(connect(c.output, pt1.input), t, systems = [pt1, c])
        sys = mtkcompile(iosys)
        prob = ODEProblem(sys, Pair[], (0.0, 100.0))
        sol = solve(prob, Rodas4())
        @test sol.retcode == Success
        @test sol[pt1.output.u]≈pt1_func.(sol.t, 1.2, 3.14) atol=1e-3

        # Test logic for a_end by constructing an integrator
        @named c = Constant(; k = 1)
        @named pt1 = TransferFunction(b = [1.2], a = [3.14, 0])
        @named iosys = System(connect(c.output, pt1.input), t, systems = [pt1, c])
        sys = mtkcompile(iosys)
        prob = ODEProblem(sys, Pair[], (0.0, 100.0))
        sol = solve(prob, Rodas4())
        @test sol.retcode == Success
        @test sol[pt1.output.u] ≈ sol.t .* (1.2 / 3.14)
        @test sol[pt1.x[1]] ≈ sol.t .* (1 / 3.14) # Test that scaling of state works properly

        # Test higher order

        function pt2_func(t, k, w, d)
            y = if d == 0
                -k * (-1 + cos(t * w))
            else
                d = complex(d)
                real(k * (1 +
                      (-cosh(sqrt(-1 + d^2) * t * w) -
                       (d * sinh(sqrt(-1 + d^2) * t * w)) / sqrt(-1 + d^2)) /
                      exp(d * t * w)))
            end
        end

        k, w, d = 1.0, 1.0, 0.5
        @named pt1 = TransferFunction(b = [w^2], a = [1, 2d * w, w^2])
        @named iosys = System(connect(c.output, pt1.input), t, systems = [pt1, c])
        sys = mtkcompile(iosys)
        prob = ODEProblem(sys, Pair[], (0.0, 100.0))
        sol = solve(prob, Rodas4())
        @test sol.retcode == Success
        @test sol[pt1.output.u]≈pt2_func.(sol.t, k, w, d) atol=1e-3

        # test zeros (high-pass version of first test)
        @named c = Constant(; k = 1)
        @named pt1 = TransferFunction(b = [1, 0], a = [1, 1])
        @named iosys = System(connect(c.output, pt1.input), t, systems = [pt1, c])
        sys = mtkcompile(iosys)
        prob = ODEProblem(sys, Pair[], (0.0, 100.0))
        sol = solve(prob, Rodas4())
        @test sol.retcode == Success
        @test sol[pt1.output.u]≈1 .- pt1_func.(sol.t, 1, 1) atol=1e-3
        @test sol[pt1.x[1]]≈pt1_func.(sol.t, 1, 1) atol=1e-3 # Test that scaling of state works properly

        # Test with no state
        @named pt1 = TransferFunction(b = [2.7], a = [pi])
        @named iosys = System(connect(c.output, pt1.input), t, systems = [pt1, c])
        sys = mtkcompile(iosys)
        prob = ODEProblem(sys, Pair[], (0.0, 100.0))
        sol = solve(prob, Rodas4())
        @test sol.retcode == Success
        @test all(==(2.7 / pi), sol[pt1.output.u])
    end
end
