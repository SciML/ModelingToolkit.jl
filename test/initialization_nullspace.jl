using ModelingToolkit, Test
using ModelingToolkit: t_nounits as t, D_nounits as D

# Underdetermined initialization: a point constrained to the unit circle, with no
# initial condition pinning where on the circle (or with what velocity) it starts.
# The initialization Jacobian is rank deficient, so the report must flag a nonempty
# null space.
@variables x(t) y(t)
@named sys = System([D(x) ~ y, 0 ~ x^2 + y^2 - 1], t)
sys = mtkcompile(sys)
prob = ODEProblem(sys, [], (0.0, 1.0); guesses = [x => 0.6, y => 0.8],
    warn_initialize_determined = false)
res = report_initialization_nullspace(prob; verbose = false)
@test res.nullity ≥ 1
@test !isempty(res.variables)
@test res.jacobian !== nothing
@test all(0 .≤ last.(res.variables) .≤ 1 + 1e-8)

# Determined initialization: an algebraic unknown uniquely fixed by a pinned state has
# a full-column-rank initialization Jacobian and no null space.
@variables u(t) w(t)
@named sys2 = System([D(u) ~ -u, 0 ~ w - exp(u)], t)
sys2 = mtkcompile(sys2)
prob2 = ODEProblem(sys2, [u => 1.0], (0.0, 1.0); guesses = [w => 1.0])
res2 = report_initialization_nullspace(prob2; verbose = false)
@test res2.nullity == 0
@test isempty(res2.variables)

# A system fully determined by its initial conditions has no initialization problem to
# analyze; the utility handles that gracefully rather than erroring.
@variables a(t) b(t)
@named sys3 = System([D(a) ~ -a, D(b) ~ -b], t)
sys3 = mtkcompile(sys3)
prob3 = ODEProblem(sys3, [a => 1.0, b => 2.0], (0.0, 1.0))
res3 = report_initialization_nullspace(prob3; verbose = false)
@test res3.nullity == 0
