using ModelingToolkitBase, OrdinaryDiffEq, StochasticDiffEq, SymbolicIndexingInterface
import Logging
using ModelingToolkitBase: t_nounits as t, D_nounits as D, ASSERTION_LOG_VARIABLE
import DiffEqNoiseProcess
import SparseConnectivityTracer as SCT
using Test

@testset "Assertion sparsity tracing" begin
    @test ModelingToolkitBase._nan_condition(true) === 0.0
    @test isnan(ModelingToolkitBase._nan_condition(false))

    @variables x y
    @named sys = System(
        [0 ~ x^2 - 1, 0 ~ y^2 - 4], [x, y], [];
        assertions = [(x * y > 0) => "Product must be positive"]
    )
    sys = complete(sys)
    residual, residual! = generate_rhs(sys; expression = Val{false})
    u = [1.0, 2.0]
    p = Float64[]
    output = similar(u)
    @test residual(u, p) == [0.0, 0.0]
    residual!(output, u, p)
    @test output == [0.0, 0.0]
    @test isnan(last(residual([-1.0, 2.0], p)))
    residual!(output, [-1.0, 2.0], p)
    @test isnan(last(output))

    detector = SCT.TracerSparsityDetector()
    # The assertion involves both variables, but contributes no derivatives.
    expected = Bool[1 0; 0 1]
    @test SCT.jacobian_sparsity(v -> residual(v, p), u, detector) == expected
    @test SCT.jacobian_sparsity((out, v) -> residual!(out, v, p), output, u, detector) ==
          expected
    @test SCT.hessian_sparsity(v -> sum(residual(v, p)), u, detector) == expected
end

@variables x(t)
@brownians a
@named inner_ode = System(D(x) ~ -sqrt(x), t; assertions = [(x > 0) => "ohno"])
@named inner_sde = System([D(x) ~ -10sqrt(x) + 0.01a], t; assertions = [(x > 0) => "ohno"])
sys_ode = mtkcompile(inner_ode)
sys_sde = mtkcompile(inner_sde)
SEED = 42

@testset "assertions are present in generated `f`" begin
    @testset "$(Problem)" for (Problem, sys, alg) in [
            (ODEProblem, sys_ode, Tsit5()), (SDEProblem, sys_sde, ImplicitEM()),
        ]
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
            (ODEProblem, sys_ode, Tsit5()), (SDEProblem, sys_sde, ImplicitEM()),
        ]
        kwargs = Problem == SDEProblem ? (; seed = SEED) : (;)
        dsys = debug_system(sys; functions = [])
        @test is_parameter(dsys, ASSERTION_LOG_VARIABLE)
        prob = Problem(dsys, [x => 0.1], (0.0, 5.0); kwargs...)
        sol = @test_logs (:error, r"ohno") match_mode = :any solve(prob, alg)
        @test !SciMLBase.successful_retcode(sol)
        prob.ps[ASSERTION_LOG_VARIABLE] = false
        sol = @test_logs min_level = Logging.Error solve(prob, alg)
        @test !SciMLBase.successful_retcode(sol)
    end
end

@testset "Hierarchical system" begin
    @testset "$(Problem)" for (ctor, Problem, inner, alg) in [
            (System, ODEProblem, inner_ode, Tsit5()),
            (System, SDEProblem, inner_sde, ImplicitEM()),
        ]
        kwargs = Problem == SDEProblem ? (; seed = SEED) : (;)
        @mtkcompile outer = ctor(Equation[], t; systems = [inner])
        dsys = debug_system(outer; functions = [])
        @test is_parameter(dsys, ASSERTION_LOG_VARIABLE)
        prob = Problem(dsys, [inner.x => 0.1], (0.0, 5.0); kwargs...)
        sol = @test_logs (:error, r"ohno") match_mode = :any solve(prob, alg)
        @test !SciMLBase.successful_retcode(sol)
        prob.ps[ASSERTION_LOG_VARIABLE] = false
        sol = @test_logs min_level = Logging.Error solve(prob, alg)
        @test !SciMLBase.successful_retcode(sol)
    end
end

@testset "Empty system doesn't error when generating assertions" begin
    @variables x(t) y(t)
    eqs = [
        0 ~ x - y
        0 ~ 2y - x
    ]

    @mtkcompile sys = System(eqs, t; assertions = [(x == 0) => "HEY!"])
    @test_nowarn generate_rhs(sys)
end
