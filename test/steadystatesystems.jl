using ModelingToolkit
using SteadyStateDiffEq
using Test
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters r
@variables x(t)
eqs = [D(x) ~ x^2 - r]
@named de = ODESystem(eqs, t)
de = complete(de)

for factor in [1e-1, 1e0, 1e10],
    u0_p in [(2.34, 2.676), (22.34, 1.632), (0.3, 15.676), (0.3, 0.006)]

    u0 = [x => factor * u0_p[1]]
    p = [r => factor * u0_p[2]]
    ss_prob = SteadyStateProblem(de, u0, p)
    sol = solve(ss_prob, SSRootfind()).u[1]
    @test abs(sol^2 - factor * u0_p[2]) < 1e-8
    ss_prob = SteadyStateProblemExpr(de, u0, p)
    sol_expr = solve(eval(ss_prob), SSRootfind()).u[1]
    @test all(x -> x == 0, sol - sol_expr)
end
