using ModelingToolkit, NonlinearSolve, SciMLBase, Test
using ModelingToolkitBase: t_nounits as t, D_nounits as D

@testset "SCC zero-state" begin
    @testset "Empty state, iip=$iip" for iip in (true, false)
        @variables x(t) y(t)
        @parameters a
        # Both unknowns compile away to observed; y depends on x.
        sys = mtkcompile(System([x ~ a, y ~ 2x + 1], t, [x, y], [a]; name = :eliminated))
        # Non-unit parameter so sol[x] == a is not a trivial 1.0 coincidence.
        prob = SCCNonlinearProblem{iip}(sys, Dict(a => 3.0))
        @test isempty(prob.u0)
        sol = solve(prob)
        @test SciMLBase.successful_retcode(sol)
        @test sol[x] == 3.0
        @test sol[y] == 7.0 # 2*3 + 1
        # remake with a different parameter; expected values derived by hand.
        sol2 = solve(remake(prob; p = [a => 5.0]))
        @test SciMLBase.successful_retcode(sol2)
        @test sol2[x] == 5.0
        @test sol2[y] == 11.0 # 2*5 + 1
    end

end
