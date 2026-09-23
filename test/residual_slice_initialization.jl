using Test, ModelingToolkit, ModelingToolkitBase, SciMLBase, OrdinaryDiffEqBDF

@testset "Residual slice initialization with algebraic unknown" begin
    @independent_variables t_residual
    @variables y(t_residual)[1:3] z(t_residual)
    D = Differential(t_residual)
    sys = complete(
        System(
            [zeros(3) ~ D(y[1:3]) + z * y[1:3], 0 ~ z - sum(y)],
            t_residual, [y, z], []; name = :residual_algebraic
        )
    )
    prob = DAEProblem(
        sys, [y => [1.0, 2.0, 3.0], D(z) => 0.0], (0.0, 1.0);
        build_initializeprob = true,
        guesses = [z => 0.0, D(y) => zeros(3)]
    )
    @test SciMLBase.successful_retcode(solve(prob.f.initializeprob))

    integ = init(prob, DFBDF())
    @test integ.u[end] ≈ 6.0
    @test integ.du[1:3] ≈ [-6.0, -12.0, -18.0]
    residual = similar(integ.u)
    prob.f(residual, integ.du, integ.u, integ.p, 0.0)
    @test residual ≈ zeros(length(residual))
end
