using ModelingToolkit, Test, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve
using LinearAlgebra
using ModelingToolkit: SciMLBase
import SciMLLogging, Logging
using ModelingToolkit: t_nounits as t, D_nounits as D

# The unknowns (by index) whose participation in mode `k` exceeds `threshold`.
participants(P, k; threshold = 0.05) = sort(findall(>(threshold), P[:, k]))

# Two decoupled linear modes with known time constants: a fast state `x` (λ = -1000) and
# a slow state `y` (λ = -1). `x` participates mostly in the fastest mode and `y` in the
# slowest.
@variables x(t) y(t)
@parameters kf = 1000.0 ks = 1.0
@named sys = System([D(x) ~ -kf * x + y, D(y) ~ -ks * y], t)
sys = mtkcompile(sys)
prob = ODEProblem(sys, [x => 1.0, y => 1.0], (0.0, 1.0))
sol = solve(prob, Rodas5P())
res = analyze_solution_jacobian(sol; verbose = false)
@test res.t == sol.t
@test length(res.eigenvalues) == length(sol.t)
λ = res.eigenvalues[1]
@test real(λ[1]) ≈ -1000 rtol = 1.0e-5
@test real(λ[end]) ≈ -1 rtol = 1.0e-5
@test res.time_constants[1] ≈ [1.0e-3, 1.0] rtol = 1.0e-5
@test all(fu -> isequal(first(only(fu)), x), res.fast_unknowns)
@test all(su -> isequal(first(only(su)), y), res.slow_unknowns)
@test all(cs -> all(==(1), cs), res.cluster_sizes)
# Participations of every mode sum to one.
@test all(P -> all(s -> isapprox(s, 1; atol = 1.0e-10), sum(P; dims = 1)), res.participation)
@test isequal(sort(res.unknowns; by = string), sort([x, y]; by = string))

# Analysis at user-selected times interpolates the solution, and the analytic Jacobian
# agrees with the default ForwardDiff Jacobian and an explicitly passed AD backend.
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

# An unstable oscillator: the modes are a complex pair with positive real part, negative
# relative damping, and participation shared by position and velocity.
@variables q(t) v(t)
@named sysu = System([D(q) ~ v, D(v) ~ -q + 0.1v], t)
sysu = mtkcompile(sysu)
solu = solve(ODEProblem(sysu, [q => 1.0, v => 0.0], (0.0, 1.0)), Rodas5P())
resu = analyze_solution_jacobian(solu; verbose = false)
@test all(l -> all(≈(0.05), real.(l)) && all(≈(sqrt(1 - 0.05^2)), abs.(imag.(l))), resu.eigenvalues)
@test all(ζ -> all(≈(-0.05), ζ), resu.damping)
@test all(ω -> all(≈(1), ω), resu.natural_frequencies)
@test all(fu -> length(fu) == 2, resu.fast_unknowns)

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
@test all(fu -> isequal(first(first(fu)), r), resd.fast_unknowns)
@test count(≈(1), resd.damping[1]) == 1
@test res.least_damped === nothing

# Conjugate pairs are listed once in the report.
report = ModelingToolkit._solution_jacobian_report(
    resd.t, resd.unknowns, resd.eigenvalues, resd.cluster_sizes, resd.time_constants,
    resd.damping, resd.participation, resd.least_damped, 0.05, 3
)
@test occursin("±", report)
@test count("ζ = 0.02", report) == 2 # the least-damped line and the mode listing
@test length(ModelingToolkit._conjugate_groups(resd.eigenvalues[1])) == 2

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
@test SciMLBase.successful_retcode(soldae)
@test length(soldae.t) > 2
resdae = analyze_solution_jacobian(soldae; verbose = false)
@test isequal(only(resdae.unknowns), a)
@test all(l -> isapprox(only(l), -101; rtol = 1.0e-5), resdae.eigenvalues)

# Problems not built from a `System` are supported; unknowns are reported by index.
plainprob = ODEProblem((u, p, t) -> [-10 * u[1], -u[2]], [1.0, 1.0], (0.0, 1.0))
resplain = analyze_solution_jacobian(solve(plainprob, Rodas5P()); verbose = false)
@test resplain.unknowns == [1, 2]
@test all(fu -> first(only(fu)) == 1, resplain.fast_unknowns)
@test all(su -> first(only(su)) == 2, resplain.slow_unknowns)

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
@test first(ModelingToolkit._rhs_jacobian(splitprob.f, [1.0], [1.0], nothing, 0.0, nothing)) ≈ fill(1.0, 1, 1) rtol = 1.0e-6

# An in-place analytic Jacobian that writes only its structurally nonzero entries is
# evaluated into a zeroed buffer.
fpart = ODEFunction(
    (du, u, p, t) -> (du[1] = -u[1]; du[2] = -2u[2]; nothing);
    jac = (J, u, p, t) -> (J[1, 1] = -1; J[2, 2] = -2; nothing), jac_prototype = zeros(2, 2)
)
respart = analyze_solution_jacobian(solve(ODEProblem(fpart, [1.0, 1.0], (0.0, 1.0)), Rodas5P()); verbose = false)
@test all(J -> J == [-1.0 0; 0 -2], respart.jacobians)

# Repeated conjugate pairs (two identical uncoupled oscillators): the coincident eigenvalues
# have no unique individual eigenvectors and form clusters involving both oscillators; their
# eigenvalues remain resolved.
λrep, Prep, κrep, crep, rrep = ModelingToolkit._modes([0 1 0 0; -4 -0.1 0 0; 0 0 0 1; 0 0 -4 -0.1])
@test length(ModelingToolkit._conjugate_groups(λrep)) == 2
@test all(==(2), crep) && all(rrep)
@test all(k -> participants(Prep, k) == [1, 2, 3, 4], 1:4)

# Two identical lightly damped oscillators coupled to a double integrator (whose defective
# eigenvalue makes the eigenbasis singular) and a better damped oscillator: the least-damped
# mode is that of the identical oscillators.
Jld = zeros(8, 8)
Jld[1:2, 1:2] = Jld[3:4, 3:4] = [0 1; -4 -0.1]
Jld[5:6, 5:6] = [0 1; -9 -3.0]
Jld[7:8, 7:8] = [0 1; 0 0.0]
Jld[8, 1] = 0.5
Jld[8, 3] = -0.5
fld = ODEFunction((u, p, t) -> Jld * u; jac = (u, p, t) -> Jld)
resld = analyze_solution_jacobian(
    solve(ODEProblem(fld, ones(8), (0.0, 0.1)), Rodas5P()); ts = [0.0], verbose = false
)
@test resld.least_damped.damping ≈ 0.025
@test sort(first.(resld.least_damped.unknowns)) == [1, 2, 3, 4]

# Nearly defective one-way coupling: the modes are individual in any units of the unknowns.
# An exactly defective pair forms a cluster in any units of the unknowns.
for c in (1.0, 1.0e-4)
    _, Pn, _, cn, _ = ModelingToolkit._modes([-1 c; 0 -1 - 1.0e-8])
    @test cn == [1, 1]
    @test sort([participants(Pn, 1), participants(Pn, 2)]) == [[1], [2]]
end
for c in (1.0, 1.0e-11)
    _, Pj, _, cj, rj = ModelingToolkit._modes([-1 c; 0 -1])
    @test cj == [2, 2] && all(rj)
    @test all(k -> participants(Pj, k) == [1, 2], 1:2)
end

# Modes are ranked by natural frequency: an undamped high-frequency oscillator is the
# fastest mode, a slow lag the slowest.
fhf = ODEFunction((u, p, t) -> [u[2], -1.0e8 * u[1], -0.01 * u[3]])
reshf = analyze_solution_jacobian(
    solve(ODEProblem(fhf, [0.0, 0.0, 1.0], (0.0, 1.0)), Rodas5P()); ts = [0.0], verbose = false
)
@test sort(first.(only(reshf.fast_unknowns))) == [1, 2]
@test first(only(only(reshf.slow_unknowns))) == 3

# Discrete-time problems are rejected.
probdisc = DiscreteProblem((u, p, t) -> 0.5 .* u, [1.0], (0.0, 3.0))
soldisc = SciMLBase.build_solution(probdisc, nothing, [0.0], [[1.0]])
@test_throws ArgumentError analyze_solution_jacobian(soldisc; verbose = false)

# The attribution does not depend on the units of the unknowns: a strong one-way coupling
# and a diagonal change of units leave the participation factors and condition numbers
# unchanged, with every mode attributed to its own unknown.
λA, PA, κA, cA, _ = ModelingToolkit._modes([-1 1.0e7; 0 -2])
@test real(λA) ≈ [-2, -1]
@test participants(PA, 1) == [2] && participants(PA, 2) == [1]
@test all(<(10), κA) && all(==(1), cA)
Dsc = Diagonal([1, 1.0e9])
A0 = [-1000 1; 1 -0.01]
λ0, P0, κ0, _, _ = ModelingToolkit._modes(A0)
λ1, P1, κ1, _, _ = ModelingToolkit._modes(Dsc * A0 / Dsc)
@test λ1 ≈ λ0
@test P1 ≈ P0 rtol = 1.0e-6
@test κ1 ≈ κ0 rtol = 1.0e-3

# A defective pair (two identical lags in series, with a strong coupling): the pair forms
# a cluster of size two whose invariant subspace involves both unknowns.
λdef, Pdef, κdef, cdef, _ = ModelingToolkit._modes([-1 1.0e6; 0 -1])
@test all(≈(-1), real(λdef))
@test all(==(2), cdef)
@test all(isfinite, Pdef) && all(k -> participants(Pdef, k) == [1, 2], 1:2)

# A double integrator in a non-modal basis: the split, ill-conditioned double zero
# eigenvalue forms a cluster of size two; the other modes are resolved individually.
Qdi = [1.0 2 0 1; 0 1 3 0; 1 0 1 2; 2 1 0 1]
λdi, Pdi, κdi, cdi, rdi = ModelingToolkit._modes(Qdi * [0 1 0 0; 0 0 0 0; 0 0 -1.0e3 0; 0 0 0 -1] / Qdi)
@test real(λdi[1]) ≈ -1000 && real(λdi[2]) ≈ -1
@test cdi == [1, 1, 2, 2] && rdi == [true, true, false, false]
@test all(isfinite, Pdi)

# Two identical defective blocks: the coincident eigenvalues form one cluster whose
# invariant subspace involves the unknowns of both blocks; every mode is clustered.
_, P22, _, c22, _ = ModelingToolkit._modes(kron(Matrix{Float64}(I, 2, 2), [-1.0 1; 0 -1]))
@test all(==(4), c22) && all(isfinite, P22)
@test all(k -> participants(P22, k) == [1, 2, 3, 4], 1:4)

# Ill-conditioned similarity transforms of well-separated eigenvalues produce several
# overlapping candidate clusters; clusters are merged entirely, so the cluster sizes are
# consistent (the number of modes with cluster size `s` is a multiple of `s`).
using Random: Xoshiro
for seed in 1:20
    Qov = ones(5, 5) .+ 1.0e-8 .* randn(Xoshiro(seed), 5, 5)
    _, Pov, _, cov, _ = ModelingToolkit._modes(Qov * Diagonal(-1.0:-1.0:-5.0) / Qov)
    @test all(s -> count(==(s), cov) % s == 0, unique(cov))
    @test all(isfinite, Pov)
end

# A defective block with exactly parallel eigenvectors does not affect the attribution of
# the isolated modes, which are attributed to their own unknowns.
Adef = diagm([0.0, 0, 0, -2, 1])
Adef[1, 2] = Adef[2, 3] = 1
λad, Pad, κad, cad, _ = ModelingToolkit._modes(Adef)
@test real(λad[1]) ≈ -2 && real(λad[2]) ≈ 1
@test participants(Pad, 1) == [4] && participants(Pad, 2) == [5]
@test cad[1:2] == [1, 1] && all(==(3), cad[3:5])
@test all(k -> participants(Pad, k) == [1, 2, 3], 3:5)

# A repeated eigenvalue with independent eigenvectors next to a defective pair: the
# repeated eigenvalue -2 is attributed to the unknowns of its eigenspace.
Arep = zeros(5, 5)
Arep[1:2, 1:2] = [-1 1; 0 -1]
Arep[3:5, 3:5] = [-2 1 0; 0 -3 0; 0 1 -2]
λr, Pr, κr, cr, _ = ModelingToolkit._modes(Arep)
for k in findall(l -> real(l) ≈ -2, λr)
    @test all(isfinite, Pr[:, k])
    @test issubset(participants(Pr, k), [3, 4, 5])
end
for k in findall(l -> real(l) ≈ -1, λr)
    @test participants(Pr, k) == [1, 2]
end

# Nilpotent matrices of order four (A^4 = 0, A^3 ≠ 0), alone and as two identical blocks:
# the split eigenvalues of each block form a cluster with finite participation.
A4 = [-3.0 6 -4 1; -1 1 0 0; 0 -1 1 0; 0 0 -1 1]
@test iszero(A4^4) && !iszero(A4^3)
_, P4, _, c4, _ = ModelingToolkit._modes(A4)
@test all(==(4), c4) && all(isfinite, P4)
_, P44, _, c44, _ = ModelingToolkit._modes(kron(Matrix{Float64}(I, 2, 2), A4))
@test all(>(1), c44) && all(isfinite, P44)

# End to end: a double integrator coupled to a fast and a slow decaying state. The slowest
# modes are the double integrator, attributed to its unknowns.
@variables xd(t) vd(t) wd(t) zd(t)
@named sysdi = System(
    [
        D(xd) ~ vd + 0.5wd, D(vd) ~ 0.3wd - 0.2zd, D(wd) ~ -1000wd + 3zd,
        D(zd) ~ -zd,
    ], t
)
sysdi = mtkcompile(sysdi)
probdi = ODEProblem(sysdi, [xd => 0.0, vd => 1.0, wd => 1.0, zd => 1.0], (0.0, 1.0))
soldi = solve(probdi, Rodas5P())
resdi = analyze_solution_jacobian(soldi; verbose = false)
@test all(fu -> isequal(first(first(fu)), wd), resdi.fast_unknowns)
@test all(su -> issetequal(Set(string.(first.(su))), Set(string.([xd, vd]))), resdi.slow_unknowns)
@test all(cs -> cs[end] == 2, resdi.cluster_sizes)

# DAE reduction: a singular algebraic block is rejected; a genuinely small algebraic row is
# accepted; the algebraic and differential equations may be scaled independently.
@test_throws ArgumentError ModelingToolkit._reduced_jacobian(
    [-1.0 1; 1 0], Diagonal([1.0, 0]), [1], 0.0, eps()
)
@test only(ModelingToolkit._reduced_jacobian([-1.0 1; 1 -2], Diagonal([1.0, 0]), [1], 0.0, eps())) ≈ -0.5
@test only(ModelingToolkit._reduced_jacobian([-1.0 1; 2.0^-40 -2.0^-39], Diagonal([1.0, 0]), [1], 0.0, eps())) ≈ -0.5
# An algebraic equation dominated by a differential unknown is not mistaken for a singular
# algebraic block.
@test only(ModelingToolkit._reduced_jacobian([-1.0 1; 1.0e7 -2], Diagonal([1.0, 0]), [1], 0.0, sqrt(eps()))) ≈ -1 + 5.0e6
for s in (1.0f0, 2.0f0^-20)
    Js = Matrix(Diagonal(Float32[-1, s, -1]))
    fms = ODEFunction(
        (du, u, p, t) -> (du .= Js * u; nothing);
        mass_matrix = Diagonal(Float32[1, s, 0]), jac = (J, u, p, t) -> (J .= Js; nothing)
    )
    solms = solve(ODEProblem(fms, Float32[1, 1.0f-3, 0], (0.0f0, 1.0f0)), Rodas5P())
    @test SciMLBase.successful_retcode(solms)
    resms = analyze_solution_jacobian(solms; verbose = false)
    @test all(l -> sort(real.(l)) ≈ [-1, 1], resms.eigenvalues)
end
# A nonsingular non-diagonal mass matrix and a purely algebraic problem.
Mnd = [1.0 1; 0 1.0e-8]
@test ModelingToolkit._reduced_jacobian(-Mnd, Mnd, [1, 2], 0.0, eps()) ≈ -I
@test size(ModelingToolkit._reduced_jacobian(ones(1, 1), zeros(1, 1), Int[], 0.0, eps())) == (0, 0)

# Solutions saved with `save_idxs` do not hold the complete state vector; for a problem
# with a symbolic system, any `save_idxs` is recorded and rejected, including a
# permutation of all unknowns.
@test_throws ArgumentError analyze_solution_jacobian(
    solve(plainprob, Rodas5P(); save_idxs = [2]); verbose = false
)
@test_throws ArgumentError analyze_solution_jacobian(
    solve(prob, Rodas5P(); save_idxs = [2, 1]); verbose = false
)

# In-place right-hand side under the default `AutoSpecialize` (function wrappers).
probiip = ODEProblem((du, u, p, t) -> (du .= -2 .* u), [1.0, 2.0], (0.0, 1.0))
resiip = analyze_solution_jacobian(solve(probiip, Rodas5P()); verbose = false)
@test all(l -> all(≈(-2), l), resiip.eigenvalues)

# Solutions without dense output (`saveat`).
probsave = ODEProblem((u, p, t) -> -u, [1.0], (0.0, 1.0))
ressave = analyze_solution_jacobian(solve(probsave, Rodas5P(); saveat = 1.0); verbose = false)
@test all(l -> only(l) ≈ -1, ressave.eigenvalues)

# A `DynamicalODEFunction` (no `f` field, partitioned state) is evaluated through its own
# callable: q' = v, v' = -4q.
probdyn = DynamicalODEProblem((dv, v, q, p, t) -> (dv .= -4 .* q), (dq, v, q, p, t) -> (dq .= v), [0.0], [1.0], (0.0, 1.0))
resdyn = analyze_solution_jacobian(solve(probdyn, Rodas5P()); verbose = false)
@test all(l -> sort(imag.(l)) ≈ [-2, 2] && all(≈(0; atol = 1.0e-10), real.(l)), resdyn.eigenvalues)

# Static arrays keep their type when evaluated with dual numbers.
SA = ModelingToolkit.StaticArrays
probsa = ODEProblem((u, p, t) -> SA.SVector(-2u[1], -3u[2]), SA.SVector(1.0, 2.0), (0.0, 1.0))
ressa = analyze_solution_jacobian(solve(probsa, Rodas5P()); verbose = false)
@test all(l -> sort(real.(l)) ≈ [-3, -2], ressa.eigenvalues)

# Complex states and Jacobians of the wrong size are rejected at the API boundary.
probc = ODEProblem((u, p, t) -> (-1 + 2im) .* u, [1.0 + 0im], (0.0, 1.0))
solc = SciMLBase.build_solution(probc, nothing, [0.0], [[1.0 + 0im]])
@test_throws ArgumentError analyze_solution_jacobian(solc; verbose = false)
fbad = ODEFunction((u, p, t) -> -u; jac = (u, p, t) -> zeros(2, 3))
probbad = ODEProblem(fbad, [1.0, 1.0], (0.0, 1.0))
solbad = SciMLBase.build_solution(probbad, nothing, [0.0], [[1.0, 1.0]])
@test_throws DimensionMismatch analyze_solution_jacobian(solbad; verbose = false)

# The report is produced through the `solution_jacobian_analysis` toggle.
@test_logs (:info, r"Solution Jacobian mode analysis") match_mode = :any analyze_solution_jacobian(sol)
@test_logs min_level = Logging.Info analyze_solution_jacobian(
    sol; verbose = MTKVerbosity(solution_jacobian_analysis = SciMLLogging.Silent)
)
