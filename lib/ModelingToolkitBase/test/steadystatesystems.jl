using ModelingToolkitBase
using NonlinearSolve
using SteadyStateDiffEq
using OrdinaryDiffEq
using Test
using ModelingToolkitBase: t_nounits as t, D_nounits as D

@parameters r
@variables x(t)
eqs = [D(x) ~ x^2 - r]
@named de = System(eqs, t)
de = complete(de)

for factor in [1.0e-1, 1.0e0, 1.0e10],
        u0_p in [(2.34, 2.676), (22.34, 1.632), (0.3, 15.676), (0.3, 0.006)]

    u0 = [x => factor * u0_p[1]]
    p = [r => factor * u0_p[2]]
    ss_prob = SteadyStateProblem(de, [u0; p])
    sol = solve(ss_prob, SSRootfind()).u[1]
    @test abs(sol^2 - factor * u0_p[2]) < 1.0e-8
    ss_prob = SteadyStateProblem(de, [u0; p])
    sol_expr = solve(ss_prob, SSRootfind()).u[1]
    @test all(x -> x == 0, sol - sol_expr)
end

@testset "Issue#4174: Algebraic equations satisfied by initialization" begin
    @parameters p d a b
    @variables A(t) X(t)
    eqs = [
        D(X) ~ p - d * X,
        a * A^2 ~ X + b,
    ]
    @mtkcompile sys = System(eqs, t)

    # ODE simulation, OK.
    u0 = [X => 0.1]
    guesses = [A => 1.0]
    ps = [p => 1.0, d => 0.5, a => 2.0, b => 16.0]
    ssprob = SteadyStateProblem(sys, [u0; ps]; guesses)
    sssol = solve(ssprob, DynamicSS(Rosenbrock23()); abstol = 1.0e-8, reltol = 1.0e-8)
    @test SciMLBase.successful_retcode(sssol)
end
