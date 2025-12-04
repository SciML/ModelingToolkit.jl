using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkitStandardLibrary.Blocks: _clamp, _dead_zone
using ModelingToolkit: inputs, unbound_inputs, bound_inputs, t_nounits as t
using OrdinaryDiffEq: ReturnCode.Success

@testset "Gain" begin
    @named c = Constant(; k = 1)
    @named gain = Gain(; k = 1)
    @named int = Integrator(; k = 1)
    @named model = System(
        [
            connect(c.output, gain.input),
            connect(gain.output, int.input)
        ],
        t, systems = [int, gain, c])

    sys = mtkcompile(model)
    prob = ODEProblem(sys, Pair[int.x => 1.0], (0.0, 1.0))
    sol = solve(prob, Rodas4())

    @test isequal(unbound_inputs(sys), [])
    @test sol.retcode == Success
    @test all(sol[c.output.u] .≈ 1)
    @test sol[int.output.u][end] ≈ 2 # expected solution after 1s
end

@testset "Feedback loop" begin
    @named c = Constant(; k = 2)
    @named gain = Gain(; k = 1)
    @named int = Integrator(; k = 1)
    @named fb = Feedback(;)
    @named model = System(
        [
            connect(c.output, fb.input1),
            connect(fb.input2, int.output),
            connect(fb.output, gain.input),
            connect(gain.output, int.input)
        ],
        t,
        systems = [int, gain, c, fb])
    sys = mtkcompile(model)

    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 100.0))

    sol = solve(prob, Rodas4())
    @test isequal(unbound_inputs(sys), [])
    @test sol.retcode == Success
    @test sol[int.output.u][end] ≈ 2 # expected solution after 1s
end

@testset "Add" begin
    @named c1 = Constant(; k = 1)
    @named c2 = Sine(; frequency = 1)
    @named add = Add(;)
    @named int = Integrator(; k = 1)
    @named model = System(
        [
            connect(c1.output, add.input1),
            connect(c2.output, add.input2),
            connect(add.output, int.input)
        ],
        t,
        systems = [int, add, c1, c2])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 1.0))
    sol = solve(prob, Rodas4())
    @test isequal(unbound_inputs(sys), [])
    @test sol.retcode == Success
    @test sol[add.output.u] ≈ 1 .+ sin.(2 * pi * sol.t)

    @testset "weights" begin
        k1 = -1
        k2 = 2
        @named add = Add(; k1 = k1, k2 = k2)
        @named model = System(
            [
                connect(c1.output, add.input1),
                connect(c2.output, add.input2),
                connect(add.output, int.input)
            ],
            t,
            systems = [int, add, c1, c2])
        sys = mtkcompile(model)
        prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 1.0))
        sol = solve(prob, Rodas4())
        @test isequal(unbound_inputs(sys), [])
        @test sol.retcode == Success
        @test sol[add.output.u] ≈ k1 .* 1 .+ k2 .* sin.(2 * pi * sol.t)
    end
end

@testset "Add3" begin
    @named c1 = Constant(; k = 1)
    @named c2 = Sine(; frequency = 1)
    @named c3 = Sine(; frequency = 2)
    @named add = Add3(;)
    @named int = Integrator(; k = 1)
    @named model = System(
        [
            connect(c1.output, add.input1),
            connect(c2.output, add.input2),
            connect(c3.output, add.input3),
            connect(add.output, int.input)
        ],
        t,
        systems = [int, add, c1, c2, c3])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 1.0))
    sol = solve(prob, Rodas4())
    @test isequal(unbound_inputs(sys), [])
    @test sol.retcode == Success
    @test sol[add.output.u] ≈ 1 .+ sin.(2 * pi * sol.t) .+ sin.(2 * pi * 2 * sol.t)

    @testset "weights" begin
        k1 = -1
        k2 = 2
        k3 = -pi
        @named add = Add3(; k1 = k1, k2 = k2, k3 = k3)
        @named model = System(
            [
                connect(c1.output, add.input1),
                connect(c2.output, add.input2),
                connect(c3.output, add.input3),
                connect(add.output, int.input)
            ],
            t,
            systems = [int, add, c1, c2, c3])
        sys = mtkcompile(model)
        prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 1.0))
        sol = solve(prob, Rodas4())
        @test isequal(unbound_inputs(sys), [])
        @test sol.retcode == Success
        @test sol[add.output.u] ≈
              k1 .* 1 .+ k2 .* sin.(2 * pi * sol.t) .+ k3 .* sin.(2 * pi * 2 * sol.t)
    end
end

@testset "Product" begin
    @named c1 = Constant(; k = 2)
    @named c2 = Sine(; frequency = 1)
    @named prod = Product(;)
    @named int = Integrator(; k = 1)
    @named model = System(
        [
            connect(c1.output, prod.input1),
            connect(c2.output, prod.input2),
            connect(prod.output, int.input)
        ],
        t,
        systems = [int, prod, c1, c2])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 1.0))
    sol = solve(prob, Rodas4())
    @test isequal(unbound_inputs(sys), [])
    @test sol.retcode == Success
    @test sol[prod.output.u] ≈ 2 * sin.(2 * pi * sol.t)
end

@testset "Power" begin
    @named c1 = Sine(; frequency = 1)
    @named c2 = Constant(; k = 2)
    @named pow = Power(;)
    @named int = Integrator(; k = 1)
    @named model = System(
        [
            connect(c1.output, pow.base),
            connect(c2.output, pow.exponent),
            connect(pow.output, int.input)
        ],
        t,
        systems = [int, pow, c1, c2])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 1.0))
    sol = solve(prob, Rodas4())
    @test isequal(unbound_inputs(sys), [])
    @test sol.retcode == Success
    @test sol[pow.output.u] ≈ sin.(2 * pi * sol.t) .^ 2
end

@testset "Modulo" begin
    @named c1 = Ramp(height = 2, duration = 1, offset = 1, start_time = 0, smooth = false)
    @named c2 = Constant(; k = 1)
    @named modl = Modulo(;)
    @named model = System(
        [
            connect(c1.output, modl.dividend),
            connect(c2.output, modl.divisor)
        ],
        t,
        systems = [modl, c1, c2])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, [], (0.0, 1.0))
    sol = solve(prob, Rodas4())
    @test isequal(unbound_inputs(sys), [])
    @test sol.retcode == Success
    @test sol[modl.remainder.u] ≈ mod.(2 * sol.t, 1)
end

@testset "UnaryMinus" begin
    @named c1 = Sine(; frequency = 1)
    @named minu = UnaryMinus(;)
    @named int = Integrator(; k = 1)
    @named model = System(
        [
            connect(c1.output, minu.input),
            connect(minu.output, int.input)
        ],
        t,
        systems = [int, minu, c1])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 1.0))
    sol = solve(prob, Rodas4())
    @test isequal(unbound_inputs(sys), [])
    @test sol.retcode == Success
    @test sol[minu.output.u] ≈ -sin.(2 * pi * sol.t)
end

@testset "Floor" begin
    @named c1 = Sine(; frequency = 1)
    @named flr = Floor(;)
    @named int = Integrator(; k = 1)
    @named model = System(
        [
            connect(c1.output, flr.input),
            connect(flr.output, int.input)
        ],
        t,
        systems = [int, flr, c1])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 1.0))
    sol = solve(prob, Rodas4())
    @test isequal(unbound_inputs(sys), [])
    @test sol.retcode == Success
    @test sol[flr.output.u] ≈ floor.(sin.(2 * pi * sol.t))
end

@testset "Ceil" begin
    @named c1 = Sine(; frequency = 1)
    @named cel = Ceil(;)
    @named int = Integrator(; k = 1)
    @named model = System(
        [
            connect(c1.output, cel.input),
            connect(cel.output, int.input)
        ],
        t,
        systems = [int, cel, c1])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 1.0))
    sol = solve(prob, Rodas4())
    @test isequal(unbound_inputs(sys), [])
    @test sol.retcode == Success
    @test sol[cel.output.u] ≈ ceil.(sin.(2 * pi * sol.t))
end

@testset "Division" begin
    @named c1 = Sine(; frequency = 1)
    @named c2 = Constant(; k = 2)
    @named div = Division(;)
    @named int = Integrator(; k = 1)
    @named model = System(
        [
            connect(c1.output, div.input1),
            connect(c2.output, div.input2),
            connect(div.output, int.input)
        ],
        t,
        systems = [int, div, c1, c2])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 1.0))
    sol = solve(prob, Rodas4())
    @test isequal(unbound_inputs(sys), [])
    @test sol.retcode == Success
    @test sol[div.output.u] ≈ sin.(2 * pi * sol.t) ./ 2
end

@testset "Abs" begin
    @named c = Sine(; frequency = 1)
    @named absb = Abs(;)
    @named int = Integrator(; k = 1)
    @named model = System(
        [
            connect(c.output, absb.input),
            connect(absb.output, int.input)
        ],
        t,
        systems = [int, absb, c])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 1.0))
    sol = solve(prob, Rodas4())
    @test isequal(unbound_inputs(sys), [])
    @test sol.retcode == Success
    @test sol[absb.output.u] ≈ abs.(sin.(2 * pi * sol.t))
end

@testset "MatrixGain" begin
    K = [1 2; 3 4]
    @named gain = MatrixGain(; K)
    K = [1, 2]
    @named gain = MatrixGain(; K)
    # TODO:
end

@testset "Sum" begin
    @named s = Sum(; input.nin = 2)
    # TODO:
end

@testset "Math" begin
    for (block, func) in [
        (Abs, abs),
        (Sign, sign),
        (Sin, sin),
        (Cos, cos),
        (Tan, tan),
        (Asin, asin),
        (Acos, acos),
        (Atan, atan),
        (Sinh, sinh),
        (Cosh, cosh),
        (Tanh, tanh),
        (Exp, exp)
    ]
        @info "testing $block"
        @named source = Sine(frequency = 1, amplitude = 0.5)
        @named b = block()
        @named int = Integrator()
        @named model = System(
            [
                connect(source.output, b.input),
                connect(b.output, int.input)
            ],
            t, systems = [int, b, source])
        sys = mtkcompile(model)
        prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 1.0))
        sol = solve(prob, Rodas4())
        @test isequal(unbound_inputs(sys), [])
        @test sol.retcode == Success
        @test sol[b.output.u] ≈ func.(sol[source.output.u])
    end

    # input must be positive
    for (block, func) in [(Sqrt, sqrt), (Log, log), (Log10, log10)]
        @info "testing $block"
        @named source = Sine(; frequency = 1, offset = 2, amplitude = 0.5)
        @named b = block()
        @named int = Integrator()
        @named model = System(
            [
                connect(source.output, b.input),
                connect(b.output, int.input)
            ],
            t, systems = [int, b, source])
        sys = mtkcompile(model)
        prob = ODEProblem(sys, Pair[int.x => 0.0, b.input.u => 2.0], (0.0, 1.0))
        sol = solve(prob, Rodas4())
        @test isequal(unbound_inputs(sys), [])
        @test sol.retcode == Success
        @test sol[b.output.u] ≈ func.(sol[source.output.u])
    end
end

@testset "Atan2" begin
    @named c1 = Sine(; frequency = 1, offset = 2)
    @named c2 = Sine(; frequency = 1, offset = 1)
    @named b = Atan2(;)
    @named int = Integrator(; k = 1)
    @named model = System(
        [
            connect(c1.output, b.input1),
            connect(c2.output, b.input2),
            connect(b.output, int.input)
        ],
        t,
        systems = [int, b, c1, c2])

    sys = mtkcompile(model)
    prob = ODEProblem(sys, Pair[int.x => 0.0, b.input1.u => 2, b.input2.u => 1], (0.0, 1.0))
    sol = solve(prob, Rodas4())

    @test isequal(unbound_inputs(sys), [])
    @test all(map(u -> u in Set([b.input1.u, b.input2.u, int.input.u]), bound_inputs(sys)))
    @test all(map(u -> u in Set([b.input1.u, b.input2.u, int.input.u]), inputs(sys)))
    @test sol.retcode == Success
    @test sol[int.input.u] ≈ atan.(sol[c1.output.u], sol[c2.output.u])
end
