using ModelingToolkit, Test
using ModelingToolkit: t_nounits as t, D_nounits as D

# Underdetermined initialization: a point constrained to the unit circle, with no
# initial condition pinning where on the circle (or with what velocity) it starts.
# The initialization Jacobian is column-rank deficient, so the analysis must flag a
# nonempty set of underdetermined unknowns.
@variables x(t) y(t)
@named sys = System([D(x) ~ y, 0 ~ x^2 + y^2 - 1], t)
sys = mtkcompile(sys)
prob = ODEProblem(sys, [], (0.0, 1.0); guesses = [x => 0.6, y => 0.8],
    warn_initialize_determined = false)
res = analyze_initialization_jacobian(prob; verbose = false)
@test res.nullity ≥ 1
@test !isempty(res.underdetermined_unknowns)
@test res.jacobian !== nothing
@test all(0 .≤ last.(res.underdetermined_unknowns) .≤ 1 + 1e-8)
@test res.rank + res.nullity == size(res.jacobian, 2)

# Overdetermined initialization with a redundant equation: two consistent initial
# conditions on a single unknown. The analysis must report a redundant equation.
@variables z(t)
@named sys2 = System([D(z) ~ -z], t; initialization_eqs = [z ~ 1.0, z^2 ~ 1.0])
sys2 = mtkcompile(sys2)
prob2 = ODEProblem(sys2, [], (0.0, 1.0); guesses = [z => 0.9],
    warn_initialize_determined = false)
res2 = analyze_initialization_jacobian(prob2; verbose = false)
@test res2.redundancy ≥ 1
@test !isempty(res2.redundant_equations)
@test all(0 .≤ last.(res2.redundant_equations) .≤ 1 + 1e-8)

# Determined initialization: an algebraic unknown uniquely fixed by a pinned state has a
# full-rank initialization Jacobian — no underdetermined unknowns and no redundant
# equations.
@variables u(t) w(t)
@named sys3 = System([D(u) ~ -u, 0 ~ w - exp(u)], t)
sys3 = mtkcompile(sys3)
prob3 = ODEProblem(sys3, [u => 1.0], (0.0, 1.0); guesses = [w => 1.0])
res3 = analyze_initialization_jacobian(prob3; verbose = false)
@test res3.nullity == 0
@test res3.redundancy == 0
@test isempty(res3.underdetermined_unknowns)
@test isempty(res3.redundant_equations)

# A system fully determined by its initial conditions has no initialization problem to
# analyze; the utility handles that gracefully rather than erroring.
@variables a(t) b(t)
@named sys4 = System([D(a) ~ -a, D(b) ~ -b], t)
sys4 = mtkcompile(sys4)
prob4 = ODEProblem(sys4, [a => 1.0, b => 2.0], (0.0, 1.0))
res4 = analyze_initialization_jacobian(prob4; verbose = false)
@test res4.nullity == 0
@test res4.redundancy == 0
