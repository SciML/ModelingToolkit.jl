using ModelingToolkit, Test, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve
using LinearAlgebra
import SciMLLogging, Logging
using ModelingToolkit: t_nounits as t, D_nounits as D

# Two decoupled linear modes with known time constants: a fast state `x` (λ = -1000) and
# a slow state `y` (λ = -1). The stiffness ratio is exactly 1000, `x` alone participates
# in the fastest mode and `y` alone in the slowest.
@variables x(t) y(t)
@parameters kf = 1000.0 ks = 1.0
@named sys = System([D(x) ~ -kf * x + y, D(y) ~ -ks * y], t)
sys = mtkcompile(sys)
prob = ODEProblem(sys, [x => 1.0, y => 1.0], (0.0, 1.0))
sol = solve(prob, Rodas5P())
res = analyze_solution_jacobian(sol; verbose = false)
@test res.t == sol.t
@test length(res.eigenvalues) == length(sol.t)
@test all(r -> isapprox(r, 1000; rtol = 1.0e-5), res.stiffness_ratio)
λ = res.eigenvalues[1]
@test real(λ[1]) ≈ -1000 rtol = 1.0e-5
@test real(λ[end]) ≈ -1 rtol = 1.0e-5
@test res.time_constants[1] ≈ [1.0e-3, 1.0] rtol = 1.0e-5
@test isequal(first(only(res.fast_unknowns)), x)
@test isequal(first(only(res.slow_unknowns)), y)
@test res.max_unstable == 0
# Participations of every mode sum to one.
@test all(P -> all(s -> isapprox(s, 1; atol = 1.0e-10), sum(P; dims = 1)), res.participation)
@test isequal(sort(res.unknowns; by = string), sort([x, y]; by = string))

# Analysis at user-selected times interpolates the solution, and the analytic Jacobian,
# finite differences and AD agree.
using ADTypes: AutoForwardDiff
res_ts = analyze_solution_jacobian(sol; ts = [0.0, 0.5], verbose = false)
@test res_ts.t == [0.0, 0.5]
res_ad = analyze_solution_jacobian(sol; ts = [0.5], verbose = false, autodiff = AutoForwardDiff())
@test res_ad.jacobians[1] ≈ res_ts.jacobians[2] rtol = 1.0e-6
probjac = ODEProblem(sys, [x => 1.0, y => 1.0], (0.0, 1.0); jac = true)
soljac = solve(probjac, Rodas5P())
res_jac = analyze_solution_jacobian(soljac; ts = [0.5], verbose = false)
@test res_jac.jacobians[1] ≈ res_ts.jacobians[2] rtol = 1.0e-6

# A state-dependent time constant: D(z) ~ -z^3 has the local eigenvalue -3z^2, which
# shrinks as `z` decays, so the analysis must track the change along the trajectory.
@variables z(t)
@named sysnl = System([D(z) ~ -z^3], t)
sysnl = mtkcompile(sysnl)
solnl = solve(ODEProblem(sysnl, [z => 2.0], (0.0, 1.0)), Rodas5P())
resnl = analyze_solution_jacobian(solnl; verbose = false)
for (u, l) in zip(solnl.u, resnl.eigenvalues)
    @test only(l) ≈ -3 * u[1]^2 rtol = 1.0e-5
end
# A single decaying mode has unit stiffness ratio.
@test all(==(1.0), resnl.stiffness_ratio)

# An unstable oscillator: the modes are a complex pair with positive real part.
@variables q(t) v(t)
@named sysu = System([D(q) ~ v, D(v) ~ -q + 0.1v], t)
sysu = mtkcompile(sysu)
solu = solve(ODEProblem(sysu, [q => 1.0, v => 0.0], (0.0, 1.0)), Rodas5P())
resu = analyze_solution_jacobian(solu; verbose = false)
@test resu.max_unstable == 2
@test all(isnan, resu.stiffness_ratio)
@test all(l -> all(≈(0.05), real.(l)) && all(≈(sqrt(1 - 0.05^2)), abs.(imag.(l))), resu.eigenvalues)
@test all(ζ -> all(≈(-0.05), ζ), resu.damping)
@test all(ω -> all(≈(1), ω), resu.natural_frequencies)

# A lightly damped oscillator coupled to a fast decaying state: the relative damping and
# natural frequency of the complex pair, and the unknowns it involves, are reported in
# `least_damped`. Real eigenvalues have relative damping one.
@variables q2(t) v2(t) r(t)
@parameters ωn = 2.0 ζ0 = 0.02
@named sysd = System(
    [D(q2) ~ v2, D(v2) ~ -ωn^2 * q2 - 2ζ0 * ωn * v2, D(r) ~ -50r + q2], t
)
sysd = mtkcompile(sysd)
probd = ODEProblem(sysd, [q2 => 1.0, v2 => 0.0, r => 0.0], (0.0, 1.0))
sold = solve(probd, Rodas5P())
resd = analyze_solution_jacobian(sold; verbose = false)
ld = resd.least_damped
@test ld !== nothing
@test ld.damping ≈ 0.02 rtol = 1.0e-5
@test ld.natural_frequency ≈ 2.0 rtol = 1.0e-5
@test isequal(sort(first.(ld.unknowns); by = string), sort([q2, v2]; by = string))
@test count(≈(1), resd.damping[1]) == 1
@test res.least_damped === nothing

# Mass-matrix DAE: the algebraic unknown is eliminated by a Schur complement. The
# algebraic equation g(w) = g(-100a) with the monotone g(s) = s^3 + s cannot be solved for
# `w` by tearing, but implies w = -100a, so the underlying ODE is D(a) ~ -101a.
@variables a(t) w(t)
@named sysdae = System([D(a) ~ -a + w, 0 ~ w^3 + w + (100a)^3 + 100a], t)
sysdae = mtkcompile(sysdae)
@test length(unknowns(sysdae)) == 2
soldae = solve(
    ODEProblem(sysdae, [a => 1.0], (0.0, 0.1); guesses = [w => -100.0]), Rodas5P()
)
resdae = analyze_solution_jacobian(soldae; verbose = false)
@test isequal(only(resdae.unknowns), a)
@test all(l -> isapprox(only(l), -101; rtol = 1.0e-5), resdae.eigenvalues)

# Problems not built from a `System` are supported; unknowns are reported by index.
plainprob = ODEProblem((u, p, t) -> [-10 * u[1], -u[2]], [1.0, 1.0], (0.0, 1.0))
resplain = analyze_solution_jacobian(solve(plainprob, Rodas5P()); verbose = false)
@test resplain.unknowns == [1, 2]
@test first(only(resplain.fast_unknowns)) == 1
@test all(r -> isapprox(r, 10; rtol = 1.0e-5), resplain.stiffness_ratio)

# Scalar problems: the right-hand side receives a number, not an array.
scalarprob = ODEProblem((u, p, t) -> -u^3, 2.0, (0.0, 1.0))
scalarsol = solve(scalarprob, Rodas5P())
resscalar = analyze_solution_jacobian(scalarsol; verbose = false)
for (u, l) in zip(scalarsol.u, resscalar.eigenvalues)
    @test only(l) ≈ -3u^2 rtol = 1.0e-5
end

# For a `SplitFunction` the analytic Jacobian covers only `f1`; the analysis must
# differentiate the complete right-hand side `f1 + f2`.
splitf = SplitFunction(
    (du, u, p, t) -> (du .= -10 .* u), (du, u, p, t) -> (du .= 11 .* u);
    jac = (J, u, p, t) -> (J .= -10)
)
splitprob = SplitODEProblem(splitf, [1.0], (0.0, 1.0))
@test ModelingToolkit._rhs_jacobian(splitprob.f, [1.0], [1.0], nothing, 0.0, nothing) ≈ fill(1.0, 1, 1) rtol = 1.0e-6

# Solutions saved with `save_idxs` do not hold the complete state vector.
plainsol_idxs = solve(plainprob, Rodas5P(); save_idxs = [2])
@test_throws ArgumentError analyze_solution_jacobian(plainsol_idxs; verbose = false)

# Conjugate pairs are listed once in the report.
report = ModelingToolkit._solution_jacobian_report(
    resd.t, resd.unknowns, resd.eigenvalues, resd.time_constants, resd.damping,
    resd.participation, resd.stiffness_ratio, resd.stiffest_index,
    [0 for _ in resd.t], resd.least_damped, 0.05, 3
)
@test occursin("±", report)
@test count("ζ = 0.02", report) == 2 # the least-damped line and the mode listing
@test length(ModelingToolkit._conjugate_groups(resd.eigenvalues[1])) == 2

# The report is produced through the `solution_jacobian_analysis` toggle.
@test_logs (:info, r"Solution Jacobian mode analysis") match_mode = :any analyze_solution_jacobian(sol)
@test_logs min_level = Logging.Info analyze_solution_jacobian(
    sol; verbose = MTKVerbosity(solution_jacobian_analysis = SciMLLogging.Silent)
)
