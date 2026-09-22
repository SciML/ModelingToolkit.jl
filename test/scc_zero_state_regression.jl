using ModelingToolkit, NonlinearSolve, SciMLBase, Test
using ModelingToolkitBase: t_nounits as t, D_nounits as D

@testset "SCC zero-state" begin
    @testset "Empty state, iip=$iip" for iip in (true, false)
        @variables x(t)
        @parameters a
        sys = mtkcompile(System([x ~ a], t, [x], [a]; name = :eliminated))
        prob = SCCNonlinearProblem{iip}(sys, Dict(a => 1.0))
        @test isempty(prob.u0)
        sol = solve(prob)
        @test SciMLBase.successful_retcode(sol)
        @test sol[x] == 1.0
    end

end
