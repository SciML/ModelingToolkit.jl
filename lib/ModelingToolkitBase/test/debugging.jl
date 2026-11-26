using ModelingToolkitBase, OrdinaryDiffEq, StochasticDiffEq, SymbolicIndexingInterface
import Logging
import SciMLLogging
using ModelingToolkitBase: t_nounits as t, D_nounits as D, ASSERTION_LOG_VARIABLE
import DiffEqNoiseProcess
import SparseConnectivityTracer as SCT
using Test

@testset "Assertion sparsity tracing" begin
    f = u -> ModelingToolkitBase._nan_condition(u[1] * u[2] > 0)
    detector = SCT.TracerSparsityDetector()
    @test SCT.jacobian_sparsity(u -> [f(u)], [1.0, 2.0], detector) == falses(1, 2)
    @test SCT.hessian_sparsity(f, [1.0, 2.0], detector) == falses(2, 2)
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

struct UnprintableError <: Exception end
struct VerboseError <: Exception end
Base.showerror(io::IO, ::VerboseError) = print(io, "λ", repeat("€", 10^6))
Base.showerror(::IO, ::UnprintableError) = error("cannot print this error")
throws_unprintable(x) = x < 0 ? throw(UnprintableError()) : x
@register_symbolic throws_unprintable(x)

@testset "Symbolic instability analysis does not throw" begin
    diagnose = ModelingToolkitBase.SciMLBase.diagnose_symbolic_instability
    @variables y(t) z(t)

    # an unused observed variable that cannot be evaluated at the last accepted state
    @mtkcompile sys = System([D(x) ~ -x, y ~ sqrt(x)], t)
    msg = diagnose(sys, [-1.0], [-1.0])
    @test occursin("evaluating sqrt(x(t)) at the last accepted state threw DomainError", msg)
    # the solver's result is returned, not the error
    prob = ODEProblem(sys, [x => -1.0], (0.0, 10.0); build_initializeprob = false)
    sol = solve(prob, Tsit5(); maxiters = 1, verbose = DEVerbosity(SciMLLogging.None()))
    @test sol.retcode == ReturnCode.MaxIters
    # and rendering the error message logs no error of its own
    sol = @test_logs min_level = Logging.Error solve(
        prob, Tsit5(); maxiters = 1, verbose = DEVerbosity()
    )
    @test sol.retcode == ReturnCode.MaxIters

    # the analysis goes on after an evaluation throws, and reports it once per equation
    @mtkcompile sys2 = System([D(x) ~ 1 / x, y ~ sqrt(x), z ~ log(1 / sqrt(x))], t)
    msg = diagnose(sys2, [-1.0e-12], [-1.0e-12])
    @test count("threw DomainError", msg) == 2
    @test occursin("division by very small value x(t)", msg)

    # an error whose message cannot be printed
    @mtkcompile usys = System([D(x) ~ -x, y ~ throws_unprintable(x)], t)
    @test occursin("threw UnprintableError", diagnose(usys, [-1.0], [-1.0]))
    # a long error message is cut
    text = ModelingToolkitBase.error_text(VerboseError())
    @test isvalid(text)
    @test endswith(text, "…")
    @test ncodeunits(text) <= ModelingToolkitBase.DIAGNOSIS_MAX_EXPRESSION_LENGTH + 4

    # assertions at a complex state
    @mtkcompile asys = System([D(x) ~ -x], t; assertions = [(x != 0) => "nonzero"])
    msg = diagnose(asys, [1.0e-12 + 1.0e-12im], [1.0e-12 + 1.0e-12im])
    @test occursin("Assertion violated", msg)
end

@testset "Symbolic instability analysis of array equations" begin
    diagnose = ModelingToolkitBase.SciMLBase.diagnose_symbolic_instability
    @variables q(t)[1:2]
    # array equations are not scalarized, so each is one term with a broadcast
    function array_system(rhs)
        @named asys = System([D(q) ~ rhs], t, collect(q), [])
        return complete(asys)
    end

    msg = diagnose(array_system(1 ./ q), [0.0, 1.0], [0.0, 1.0])
    @test occursin("division by very small values of q(t) ([1]: 0)", msg)
    @test isempty(diagnose(array_system(1 ./ q), [1.0, 1.0], [1.0, 1.0]))
    msg = diagnose(array_system(sqrt.(q)), [0.0, 1.0], [0.0, 1.0])
    @test occursin("sqrt of q(t) ([1]: 0)", msg)
    msg = diagnose(array_system(log.(q)), [1.0, 0.0], [1.0, 0.0])
    @test occursin("log of q(t) ([2]: 0)", msg)
    msg = diagnose(array_system(q .^ 2), [2.0, 2.0], [2.0, 2.0])
    @test occursin("([1]: base ≈ 2, power 2; [2]: base ≈ 2, power 2)", msg)
    msg = diagnose(array_system(q .^ -1), [0.0, 1.0], [0.0, 1.0])
    @test occursin("raised to negative powers with base going to 0 ([1]: base ≈ 0, power -1)", msg)
    msg = diagnose(array_system(inv.(q)), [1.0, 0.0], [1.0, 0.0])
    @test occursin("division by very small values of q(t) ([2]: 0)", msg)
    # only the denominator needs a value
    @parameters c
    msg = diagnose(array_system(c ./ q), [0.0, 1.0], [0.0, 1.0])
    @test occursin("division by very small values of q(t) ([1]: 0)", msg)

    # an array unknown takes as many entries of the state as it has elements
    @variables x(t)
    @named xsys = System([D(q) ~ 1 ./ q, D(x) ~ 1 / x], t)
    xsys = complete(xsys)
    @test isequal(unknowns(xsys), [q, x])
    msg = diagnose(xsys, [1.0, 0.0, 1.0], [1.0, 0.0, 1.0])
    @test occursin("division by very small values of q(t) ([2]: 0)", msg)
    @test !occursin("x(t) ≈", msg)
    # an array unknown with offset axes
    @variables k(t)[0:1]
    @named ksys = System([D(k) ~ -k, D(x) ~ 1 / x], t)
    msg = diagnose(complete(ksys), [1.0, 0.0, 1.0], [1.0, 0.0, 1.0])
    @test !occursin("x(t) ≈", msg)

    # an array term with more elements than the budget left is not checked, and stops
    # the analysis
    ctx = ModelingToolkitBase.SingularityAnalysis(xsys, Dict{ModelingToolkitBase.SymbolicT, ModelingToolkitBase.SymbolicT}(), String[])
    ctx.budget[] = 1000
    singular = (b, e) -> e > 0 && abs(b) > 1
    vals = Any[fill(2.0, 1001), ones(1, 1001)]
    @test ModelingToolkitBase.failing_elements(singular, string, vals, ctx) === nothing
    @test ctx.budget[] == 0

    # a term lists at most `DIAGNOSIS_MAX_ELEMENTS` of its elements
    n = ModelingToolkitBase.DIAGNOSIS_MAX_ELEMENTS + 3
    @variables p(t)[1:n]
    @named psys = System([D(p) ~ 1 ./ p], t, collect(p), [])
    msg = diagnose(complete(psys), zeros(n), zeros(n))
    @test count("division by very small values", msg) == 1
    @test occursin("; and 3 more)", msg)

    # array-valued observed variables, after scalarization
    @variables w(t)[1:2]
    @mtkcompile osys = System([D(q) ~ 1 ./ w, w ~ q .- 1], t)
    msg = diagnose(osys, [1.0, 2.0], [1.0, 2.0])
    @test count("division by very small value", msg) == 1
    @test occursin("(w(t))[1]", msg)

    # an array-valued observed variable used by an unscalarized equation
    @named wsys = System([D(q) ~ 1 ./ w], t, collect(q), []; observed = [w ~ q .- 1])
    wsys = complete(wsys)
    msg = diagnose(wsys, [1.0, 2.0], [1.0, 2.0])
    @test occursin("division by very small values of w(t) ([1]: 0)", msg)
    @test isempty(diagnose(wsys, [3.0, 2.0], [3.0, 2.0]))

    # the analysis does not grow the system's IR
    ir = ModelingToolkitBase.get_irstructure(osys)
    nir = length(ir)
    for k in 1:5
        diagnose(osys, [k, 2.0], [k, 2.0])
    end
    @test length(ir) == nir
end

@testset "Abbreviation keeps whole characters" begin
    @test isvalid(ModelingToolkitBase.abbreviated(repeat("λ", 99) * "€" * repeat("λ", 20)))
end
