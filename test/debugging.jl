using ModelingToolkit, OrdinaryDiffEq, StochasticDiffEq, SymbolicIndexingInterface
using ModelingToolkit: t_nounits as t, D_nounits as D, ASSERTION_LOG_VARIABLE

@variables x(t)
@brownian a
@named inner_ode = ODESystem(D(x) ~ -sqrt(x), t; assertions = [(x > 0) => "ohno"])
@named inner_sde = System([D(x) ~ -sqrt(x) + a], t; assertions = [(x > 0) => "ohno"])
sys_ode = structural_simplify(inner_ode)
sys_sde = structural_simplify(inner_sde)

@testset "`debug_system` adds assertions" begin
    @testset "$(typeof(sys))" for (Problem, sys, alg) in [
        (ODEProblem, sys_ode, Tsit5()), (SDEProblem, sys_sde, ImplicitEM())]
        dsys = debug_system(sys; functions = [])
        @test is_parameter(dsys, ASSERTION_LOG_VARIABLE)
        prob = Problem(dsys, [x => 1.0], (0.0, 5.0))
        sol = solve(prob, alg)
        @test !SciMLBase.successful_retcode(sol)
        prob.ps[ASSERTION_LOG_VARIABLE] = true
        sol = @test_logs (:error, r"ohno") match_mode=:any solve(prob, alg)
        @test !SciMLBase.successful_retcode(sol)
    end
end

@testset "Hierarchical system" begin
    @testset "$(typeof(inner))" for (ctor, Problem, inner, alg) in [
        (ODESystem, ODEProblem, inner_ode, Tsit5()),
        (System, SDEProblem, inner_sde, ImplicitEM())]
        @mtkbuild outer = ctor(Equation[], t; systems = [inner])
        dsys = debug_system(outer; functions = [])
        @test is_parameter(dsys, ASSERTION_LOG_VARIABLE)
        prob = Problem(dsys, [inner.x => 1.0], (0.0, 5.0))
        sol = solve(prob, alg)
        @test !SciMLBase.successful_retcode(sol)
        prob.ps[ASSERTION_LOG_VARIABLE] = true
        sol = @test_logs (:error, r"ohno") match_mode=:any solve(prob, alg)
        @test !SciMLBase.successful_retcode(sol)
    end
end
