using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit
using OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D

@testset "Array Guesses" begin
    for (block, guess) in [
        (RealInputArray(; nin = 3, name = :a), zeros(3)),
        (RealOutputArray(; nout = 3, name = :a), zeros(3))
    ]
        guesses = ModelingToolkit.guesses(block)
        @test guesses[@nonamespace block.u] == guess
    end
end

@testset "Scalarized Guesses" begin
    for (block, guess) in [
        (RealInput(; name = :a), 0.0),
        (RealInput(; nin = 3, name = :a), zeros(3)),
        (RealOutput(; name = :a), 0.0),
        (RealOutput(; nout = 3, name = :a), zeros(3))
    ]
        guesses = ModelingToolkit.guesses(block)
        @test guesses[@nonamespace block.u[1]] == guess[1]
    end
end

@testset "SISO Check" begin
    k, w, d = 1.0, 1.0, 0.5
    @named c = Constant(; k = 1)
    @named so = SecondOrder(; k = k, w = w, d = d, xd = 1)
    @named iosys = System(connect(c.output, so.input), t, systems = [so, c])
    sys = mtkcompile(iosys)

    initsys = ModelingToolkit.generate_initializesystem(sys)
    initsys = mtkcompile(initsys)
    initprob = NonlinearProblem(initsys, [t => 0])
    initsol = solve(initprob)

    @test initsol[sys.so.xd] == 1.0
    @test initsol[sys.so.u] == 1.0
end

@test_deprecated RealInput(; name = :a, u_start = 1.0)
@test_deprecated RealInput(; name = :a, nin = 2, u_start = ones(2))
@test_deprecated RealOutput(; name = :a, u_start = 1.0)
@test_deprecated RealOutput(; name = :a, nout = 2, u_start = ones(2))
@test_deprecated RealInputArray(; name = :a, nin = 2, u_start = ones(2))
@test_deprecated RealOutputArray(; name = :a, nout = 2, u_start = ones(2))
