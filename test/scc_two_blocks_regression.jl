using ModelingToolkit, NonlinearSolve, SciMLBase, Test
using ModelingToolkitBase: t_nounits as t, D_nounits as D

@testset "SCC two-blocks" begin
    @testset "Independent variable, iip=$iip, combine=$combine" for iip in (true, false),
            combine in (true, false)
        @variables y(t) z(t)
        @parameters a b
        sys = mtkcompile(
            System([0 ~ y^2 - a, 0 ~ z^2 - b], t, [y, z], [a, b]; name = :algebraic)
        )
        prob = SCCNonlinearProblem{iip}(
            sys, Dict(a => 9.0, b => 16.0); combine_sccs = combine,
            guesses = Dict(y => 2.5, z => 3.5)
        )
        sol = solve(prob; abstol = 1.0e-12, reltol = 1.0e-12)
        @test SciMLBase.successful_retcode(sol)
        @test sol[y] ≈ 3.0
        @test sol[z] ≈ 4.0
    end

end
