using ModelingToolkit, NonlinearSolve, SciMLBase, Test
using ModelingToolkitBase: t_nounits as t, SystemCompatibilityError

@testset "SCC two-blocks" begin
    @variables y(t) z(t)
    @parameters a b
    timedep = mtkcompile(
        System([0 ~ y^2 - a, 0 ~ z^2 - b], t, [y, z], [a, b]; name = :algebraic)
    )
    op = Dict(a => 9.0, b => 16.0)
    guesses = Dict(y => 2.5, z => 3.5)

    @testset "Time-dependent system fails check_compatibility" begin
        @test_throws SystemCompatibilityError SCCNonlinearProblem(
            timedep, op; combine_sccs = false, guesses
        )
        @test_throws SystemCompatibilityError SCCNonlinearProblem(
            timedep, op; combine_sccs = true, guesses
        )
    end

    @testset "After NonlinearSystem conversion, iip=$iip, combine=$combine" for iip in (
                true, false,
            ),
            combine in (true, false)
        sys = mtkcompile(NonlinearSystem(timedep))
        prob = SCCNonlinearProblem{iip}(
            sys, op; combine_sccs = combine, guesses
        )
        sol = solve(prob; abstol = 1.0e-12, reltol = 1.0e-12)
        @test SciMLBase.successful_retcode(sol)
        @test sol[y] ≈ 3.0
        @test sol[z] ≈ 4.0
    end
end
