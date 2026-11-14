using ModelingToolkitBase, OrdinaryDiffEq, StochasticDiffEq, SymbolicIndexingInterface
import Logging
using ModelingToolkitBase: t_nounits as t, D_nounits as D, ASSERTION_LOG_VARIABLE
import DiffEqNoiseProcess
using Test

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

@testset "Symbolic instability analysis stays bounded" begin
    diagnose = ModelingToolkitBase.SciMLBase.diagnose_symbolic_instability

    # a singularity through an observed variable is found, and the message names the
    # (small) equation it appears in
    @variables y(t)
    @mtkcompile sys = System([D(x) ~ 1 / y, y ~ x - 1], t)
    msg = diagnose(sys, [1.0], [1.0])
    @test occursin("division by very small value", msg)
    @test occursin("y(t)", msg)
    @test isempty(diagnose(sys, [3.0], [3.0]))

    # subexpressions with parameters cannot be evaluated without their values, and are
    # skipped rather than substituted into
    @parameters p = 1.0
    @mtkcompile psys = System([D(x) ~ x / (p - x)], t)
    @test isempty(diagnose(psys, [1.0], [1.0]))

    # a deep chain of observed equations, each using the previous one three times: its
    # `full_equations` print exponentially large
    function chain(N)
        @variables z(t) = 1.0 w(t)[1:N]
        @parameters c[1:N] = collect(range(0.5, 1.5; length = N))
        eqs = [
            D(z) ~ -z + w[N]; w[1] ~ z / (1 + c[1] * z^2);
            [w[i] ~ (w[i - 1] + c[i]) / (1 + c[i] * w[i - 1]^2) + sqrt(c[i] + w[i - 1]^2) for i in 2:N]
        ]
        return mtkcompile(System(eqs, t; name = :chain))
    end
    csys = chain(200)
    stats = @timed diagnose(csys, [3.0], [3.0])
    @test stats.time < 5
    @test count("raised to power", stats.value) == 1
    long = ModelingToolkitBase.abbreviated(only(full_equations(chain(20))))
    @test length(long) <= ModelingToolkitBase.DIAGNOSIS_MAX_EXPRESSION_LENGTH + 2
    @test endswith(long, "…")

    # at most `DIAGNOSIS_MAX_FINDINGS` findings are listed
    n = ModelingToolkitBase.DIAGNOSIS_MAX_FINDINGS + 10
    @variables q(t)[1:n]
    @mtkcompile qsys = System([D(q[i]) ~ -q[i]^2 for i in 1:n], t)
    msg = diagnose(qsys, fill(2.0, n), fill(2.0, n))
    @test count("raised to power", msg) == ModelingToolkitBase.DIAGNOSIS_MAX_FINDINGS
    @test occursin("(and 10 more of these)", msg)
end
