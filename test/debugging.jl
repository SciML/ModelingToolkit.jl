using ModelingToolkit, OrdinaryDiffEq, StochasticDiffEq, SymbolicIndexingInterface
import Logging
using ModelingToolkit: t_nounits as t, D_nounits as D, ASSERTION_LOG_VARIABLE

@variables x(t)
@brownians a
@named inner_ode = System(D(x) ~ -sqrt(x), t; assertions = [(x > 0) => "ohno"])
@named inner_sde = System([D(x) ~ -10sqrt(x) + 0.01a], t; assertions = [(x > 0) => "ohno"])
sys_ode = mtkcompile(inner_ode)
sys_sde = mtkcompile(inner_sde)
SEED = 42

@testset "assertions are present in generated `f`" begin
    @testset "$(Problem)" for (Problem, sys, alg) in [
        (ODEProblem, sys_ode, Tsit5()), (SDEProblem, sys_sde, ImplicitEM())]
        kwargs = Problem == SDEProblem ? (; seed = SEED) : (;)
        @test !is_parameter(sys, ASSERTION_LOG_VARIABLE)
        prob = Problem(sys, [x => 0.1], (0.0, 5.0); kwargs...)
        sol = solve(prob, alg)
        @test !SciMLBase.successful_retcode(sol)
        @test isnan(prob.f.f([0.0], prob.p, sol.t[end])[1])
    end
end

@testset "`debug_system` adds logging" begin
    @testset "$(Problem)" for (Problem, sys, alg) in [
        (ODEProblem, sys_ode, Tsit5()), (SDEProblem, sys_sde, ImplicitEM())]
        kwargs = Problem == SDEProblem ? (; seed = SEED) : (;)
        dsys = debug_system(sys; functions = [])
        @test is_parameter(dsys, ASSERTION_LOG_VARIABLE)
        prob = Problem(dsys, [x => 0.1], (0.0, 5.0); kwargs...)
        sol = @test_logs (:error, r"ohno") match_mode=:any solve(prob, alg)
        @test !SciMLBase.successful_retcode(sol)
        prob.ps[ASSERTION_LOG_VARIABLE] = false
        sol = @test_logs min_level=Logging.Error solve(prob, alg)
        @test !SciMLBase.successful_retcode(sol)
    end
end

@testset "Hierarchical system" begin
    @testset "$(Problem)" for (ctor, Problem, inner, alg) in [
        (System, ODEProblem, inner_ode, Tsit5()),
        (System, SDEProblem, inner_sde, ImplicitEM())]
        kwargs = Problem == SDEProblem ? (; seed = SEED) : (;)
        @mtkcompile outer = ctor(Equation[], t; systems = [inner])
        dsys = debug_system(outer; functions = [])
        @test is_parameter(dsys, ASSERTION_LOG_VARIABLE)
        prob = Problem(dsys, [inner.x => 0.1], (0.0, 5.0); kwargs...)
        sol = @test_logs (:error, r"ohno") match_mode=:any solve(prob, alg)
        @test !SciMLBase.successful_retcode(sol)
        prob.ps[ASSERTION_LOG_VARIABLE] = false
        sol = @test_logs min_level=Logging.Error solve(prob, alg)
        @test !SciMLBase.successful_retcode(sol)
    end
end
