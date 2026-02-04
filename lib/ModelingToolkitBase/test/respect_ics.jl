using ModelingToolkit, NonlinearSolve
using ModelingToolkit: t_nounits as t, D_nounits as D
using Test

@testset "ensure ICs are respected" begin
    @parameters g
    @variables x(t) y(t) [state_priority = 10] λ(t)
    eqs = [D(D(x)) ~ λ * x
           D(D(y)) ~ λ * y - g
           x^2 + y^2 ~ 1]
    @mtkcompile pend = System(eqs, t)

    iprob = ModelingToolkit.InitializationProblem(
        pend, 0.0,
        [x => 1.0, D(y) => 0.0, g => 1],
        guesses = [λ => 1, y => 0.0])
    isol = solve(iprob)

    # previous behavior was `==` for these
    @test isol.retcode != ReturnCode.Stalled
    @test isol[isol.prob.f.sys.x] != 0.0
    @test iprob.ps[Initial(x)] != 0.0

    # current behavior
    @test isol.retcode == ReturnCode.Success
    @test isol[isol.prob.f.sys.x] == 1.0
    @test iprob.ps[Initial(x)] == 1.0
end
