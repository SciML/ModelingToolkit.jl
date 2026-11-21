using ModelingToolkit, Test, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve
using LinearAlgebra
using ModelingToolkit: SciMLBase
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
@test SciMLBase.successful_retcode(soldae)
@test length(soldae.t) > 2
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
@test first(ModelingToolkit._rhs_jacobian(splitprob.f, [1.0], [1.0], nothing, 0.0, nothing)) ≈ fill(1.0, 1, 1) rtol = 1.0e-6

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

# Classification of the eigenvalues of a matrix with entries accurate to `ϵ`.
function classify(J)
    λ, P, κ, V = ModelingToolkit._modes(J, eps())
    r = ModelingToolkit._eigenvalue_error_bounds(λ, κ, V, eps(), norm(J))
    zm = abs.(λ) .<= r
    σ = [zm[k] || abs(real(λ[k])) <= r[k] ? 0.0 : real(λ[k]) for k in eachindex(λ)]
    return (; λ, P, κ, zm, σ)
end

# Repeated conjugate pairs (two identical uncoupled oscillators) form two oscillatory
# modes, each with the participation of its own position and velocity.
crep = classify([0 1 0 0; -4 -0.1 0 0; 0 0 0 1; 0 0 -4 -0.1])
λrep, Prep = crep.λ, crep.P
@test length(ModelingToolkit._conjugate_groups(λrep)) == 2
@test all(k -> λrep[2k] == conj(λrep[2k - 1]) && imag(λrep[2k - 1]) > 0, 1:2)
@test all(k -> Prep[:, 2k] ≈ Prep[:, 2k - 1], 1:2)
@test all(k -> count(>(0.1), Prep[:, 2k]) == 2, 1:2)

# A fast real mode does not hide a slower, lightly damped oscillatory mode.
cstiff = classify([-1.0e6 0 0; 0 -0.1 10; 0 -10 -0.1])
@test !any(cstiff.zm)
ldstiff = ModelingToolkit._least_damped_mode(
    [0.0], [:a, :b, :c], [cstiff.λ], [BitVector(cstiff.zm)],
    [-cstiff.σ ./ abs.(cstiff.λ)], [cstiff.P], 0.05
)
@test ldstiff !== nothing
@test ldstiff.natural_frequency ≈ hypot(0.1, 10)
@test ldstiff.damping ≈ 0.1 / hypot(0.1, 10)
@test ModelingToolkit._stiffness_ratio(cstiff.σ) ≈ 1.0e7

# Well-conditioned small eigenvalues are resolved, however small relative to the fastest
# mode: the exact zero is the only zero eigenvalue.
cdiag = classify(Matrix(Diagonal([-1, 0, -1.0e-8, 1.0e-8])))
@test count(cdiag.zm) == 1 && iszero(only(cdiag.λ[cdiag.zm]))
@test count(>(0), cdiag.σ) == 1
@test ModelingToolkit._stiffness_ratio(cdiag.σ) ≈ 1.0e8

# A double integrator in a non-modal basis: the double zero eigenvalue is perturbed by
# about √ϵ to a pair of ill-conditioned small eigenvalues. Their condition numbers identify
# them as two zero eigenvalues rather than a slow decaying and a spurious unstable mode;
# their participation is not resolved.
Qdi = [1.0 2 0 1; 0 1 3 0; 1 0 1 2; 2 1 0 1]
cdi = classify(Qdi * [0 1 0 0; 0 0 0 0; 0 0 -1.0e3 0; 0 0 0 -1] / Qdi)
@test count(cdi.zm) == 2
@test count(>(0), cdi.σ) == 0
@test ModelingToolkit._stiffness_ratio(cdi.σ) ≈ 1.0e3 rtol = 1.0e-6
@test all(isnan, cdi.P[:, cdi.zm])

# An exact integrator next to a double integrator in a non-modal basis.
Qmix = [1.0 2 0 1 0; 0 1 3 0 1; 1 0 1 2 0; 2 1 0 1 1; 0 1 1 0 2]
cmix = classify(Qmix * [0 1 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 -1 0; 0 0 0 0 -10] / Qmix)
@test count(cmix.zm) == 3

# Nilpotent matrices of order four (A^4 = 0, A^3 ≠ 0) and three in non-modal bases: the
# eigenvalues are split by about ϵ^(1/4) and ϵ^(1/3) and must all be zero eigenvalues, not
# resolved growing modes.
A4 = [-3.0 6 -4 1; -1 1 0 0; 0 -1 1 0; 0 0 -1 1]
@test iszero(A4^4) && !iszero(A4^3)
c4 = classify(A4)
@test all(c4.zm)
@test count(>(0), c4.σ) == 0
# Two identical nilpotent blocks of order four: the coincident split eigenvalues of the two
# blocks do not form defective pairs.
c44 = classify(kron(Matrix{Float64}(I, 2, 2), A4))
@test all(c44.zm)
@test count(>(0), c44.σ) == 0
Q3 = [1.0 2 0; 0 1 3; 1 0 1]
c3 = classify(Q3 * [0 1 0; 0 0 1; 0 0 0] / Q3)
@test all(c3.zm)

# A nearly defective pair near -1: the first-order bound reaches the other eigenvalue and is
# limited to the bound of a double eigenvalue, so the pair is resolved as decaying.
cnear = classify([-1.0 1; 2.0e-28 -1])
@test !any(cnear.zm)
@test all(≈(-1), cnear.σ)

# An exactly defective pair next to well-separated modes: the first-order bound of the pair
# covers the whole spectrum, but the cluster consists of the pair only, so all modes are
# resolved as decaying.
Jsep = Matrix(Diagonal([-1.0, -1, -2, -3, -4, -5]))
Jsep[1, 2] = 1
csep = classify(Jsep)
@test !any(csep.zm)
@test sort(csep.σ) ≈ [-5, -4, -3, -2, -1, -1]

# A defective eigenvalue away from the origin is neither a zero eigenvalue nor resolved
# into participations of individual unknowns.
cdef = classify([-1.0 1; 0 -1])
@test !any(cdef.zm)
@test all(isnan, cdef.P)
@test all(≈(-1), cdef.σ)

# End to end: a double integrator coupled to a fast and a slow decaying state.
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
@test all(m -> count(m) == 2, resdi.zero_eigenvalues)
@test resdi.max_unstable == 0
@test all(r -> isapprox(r, 1000; rtol = 1.0e-6), resdi.stiffness_ratio)
@test all(d -> all(isnan, d[resdi.zero_eigenvalues[1]]), resdi.damping)
@test all(τ -> count(isinf, τ) == 2, resdi.time_constants)

# Diagonal Jacobian with slow modes of both signs (analytic Jacobian): one growing mode,
# a stiffness ratio of 1e8 and finite time constants.
fdiag = ODEFunction((u, p, t) -> [-u[1], -1.0e-8 * u[2], 1.0e-8 * u[3]]; jac = (u, p, t) -> Matrix(Diagonal([-1, -1.0e-8, 1.0e-8])))
soldiag = solve(ODEProblem(fdiag, [1.0, 1.0, 1.0], (0.0, 1.0)), Rodas5P())
resdiag = analyze_solution_jacobian(soldiag; verbose = false)
@test resdiag.max_unstable == 1
@test all(r -> r ≈ 1.0e8, resdiag.stiffness_ratio)
@test all(τ -> sort(τ) ≈ [1, 1.0e8, 1.0e8], resdiag.time_constants)

# An exact zero eigenvalue at an equilibrium is not reported as a growing mode: the default
# Jacobian (ForwardDiff) is exact there, also for single-precision states.
for T in (Float64, Float32)
    probeq = ODEProblem((u, p, t) -> [-u[1], u[2]^3], T[1, 0], (T(0), T(1)))
    reseq = analyze_solution_jacobian(solve(probeq, Rodas5P()); verbose = false)
    @test reseq.max_unstable == 0
    @test all(m -> count(m) == 1, reseq.zero_eigenvalues)
end

# A singular algebraic block is rejected, also when it is singular only to the accuracy
# of the Jacobian.
@test_throws ArgumentError ModelingToolkit._reduced_jacobian(
    [-1.0 1; 1 0], Diagonal([1.0, 0]), [1], 0.0, eps()
)
@test_throws ArgumentError ModelingToolkit._reduced_jacobian(
    [-1.0 1; 1 1.0e-12], Diagonal([1.0, 0]), [1], 0.0, sqrt(eps())
)
Jred, _ = ModelingToolkit._reduced_jacobian([-1.0 1; 1 -2], Diagonal([1.0, 0]), [1], 0.0, eps())
@test only(Jred) ≈ -0.5

# Cancellation in the Schur complement: both residuals are u[1] + 49u[2], so the reduced
# Jacobian is 1 - 49 / 49 = 0 exactly, and is computed as a rounding error. The error scale
# of the terms entering the complement classifies it as a zero eigenvalue rather than a
# growing mode, with the analytic Jacobian and with ForwardDiff.
fcancel(du, u, p, t) = (du[1] = u[1] + 49u[2]; du[2] = u[1] + 49u[2]; nothing)
for jac in (nothing, (J, u, p, t) -> (J .= [1.0 49; 1 49]))
    fc = ODEFunction(fcancel; mass_matrix = Diagonal([1.0, 0.0]), jac)
    solcancel = solve(ODEProblem(fc, [49.0, -1.0], (0.0, 1.0)), Rodas5P())
    @test SciMLBase.successful_retcode(solcancel)
    rescancel = analyze_solution_jacobian(solcancel; verbose = false)
    @test rescancel.max_unstable == 0
    @test all(only, rescancel.zero_eigenvalues)
end

# Scaling an algebraic equation does not change the reduced Jacobian or the rank test:
# x' = -x + y, 0 = s (x - 2y) has the reduced rate -0.5 for any s ≠ 0, also in Float32
# with an analytic Jacobian.
for s in (1.0f0, 2.0f0^-20)
    fs(du, u, p, t) = (du[1] = -u[1] + u[2]; du[2] = s * (u[1] - 2u[2]); nothing)
    js(J, u, p, t) = (J .= [-1 1; s -2s]; nothing)
    ffs = ODEFunction(fs; mass_matrix = Diagonal(Float32[1, 0]), jac = js)
    sols = solve(ODEProblem(ffs, Float32[1, 0.5], (0.0f0, 1.0f0)), Rodas5P())
    @test SciMLBase.successful_retcode(sols)
    ress = analyze_solution_jacobian(sols; verbose = false)
    @test all(l -> only(l) ≈ -0.5, ress.eigenvalues)
end

# A nonsingular non-diagonal mass matrix: `M⁻¹ J = -I` exactly, decaying modes.
Mnd = [1.0 1; 0 1.0e-8]
Jnd, Snd = ModelingToolkit._reduced_jacobian(-Mnd, Mnd, [1, 2], 0.0, eps())
@test Jnd ≈ -I
λnd, _, κnd, Vnd = ModelingToolkit._modes(Jnd, eps())
@test !any(abs.(λnd) .<= ModelingToolkit._eigenvalue_error_bounds(λnd, κnd, Vnd, eps(), Snd))

# A purely algebraic mass-matrix problem reduces to an empty Jacobian.
Jempty, _ = ModelingToolkit._reduced_jacobian(ones(1, 1), zeros(1, 1), Int[], 0.0, eps())
@test size(Jempty) == (0, 0)

# A finite-difference Jacobian is accurate only relative to the overall scale: an algebraic
# row of the size of rounding errors (0 = y^3 - 8x^3 at the origin) is not normalized into a
# valid pivot. A genuinely small algebraic row of an exact Jacobian is accepted.
@test_throws ArgumentError ModelingToolkit._reduced_jacobian(
    [-3.0 1; -8eps() eps()], Diagonal([1.0, 0]), [1], 0.0, sqrt(eps()), false
)
Jsmall, _ = ModelingToolkit._reduced_jacobian(
    [-1.0 1; 2.0^-40 -2.0^-39], Diagonal([1.0, 0]), [1], 0.0, eps(), true
)
@test only(Jsmall) ≈ -0.5

# A defective block does not hide isolated modes: the eigenbasis is singular, but the
# growing mode and the decaying mode with time constant 0.5 are resolved.
Adef = diagm([0.0, 0, 0, -2, 1])
Adef[1, 2] = Adef[2, 3] = 1
fdef = ODEFunction((u, p, t) -> Adef * u; jac = (u, p, t) -> Adef)
soldef3 = solve(ODEProblem(fdef, [0.0, 0, 0, 1, 1.0e-3], (0.0, 1.0)), Rodas5P())
resdef3 = analyze_solution_jacobian(soldef3; verbose = false)
@test resdef3.max_unstable == 1
@test all(m -> count(m) == 3, resdef3.zero_eigenvalues)
@test all(τ -> any(≈(0.5), τ) && any(≈(1), τ), resdef3.time_constants)

# Scaling a differential equation does not affect the resolution of the other modes:
# M = diag(1, s, 0), J = diag(-1, s, -1) has the modes -1 and +1 for any s ≠ 0.
for s in (1.0f0, 2.0f0^-20)
    Js = Matrix(Diagonal(Float32[-1, s, -1]))
    fms = ODEFunction(
        (du, u, p, t) -> (du .= Js * u; nothing);
        mass_matrix = Diagonal(Float32[1, s, 0]), jac = (J, u, p, t) -> (J .= Js; nothing)
    )
    solms = solve(ODEProblem(fms, Float32[1, 1.0f-3, 0], (0.0f0, 1.0f0)), Rodas5P())
    @test SciMLBase.successful_retcode(solms)
    resms = analyze_solution_jacobian(solms; verbose = false)
    @test resms.max_unstable == 1
    @test all(l -> sort(real.(l)) ≈ [-1, 1], resms.eigenvalues)
    @test all(τ -> all(isfinite, τ), resms.time_constants)
end

# A defective Jacobian: no participation and no fast/slow unknowns are reported.
soldef = solve(ODEProblem((u, p, t) -> [-u[1] + u[2], -u[2]], [1.0, 1.0], (0.0, 1.0)), Rodas5P())
resdef = analyze_solution_jacobian(soldef; verbose = false)
@test all(P -> all(isnan, P), resdef.participation)
@test isempty(resdef.fast_unknowns) && isempty(resdef.slow_unknowns)

# For a problem with a symbolic system, any `save_idxs` is recorded in the solution and
# rejected, including a permutation of all unknowns.
@test_throws ArgumentError analyze_solution_jacobian(
    solve(prob, Rodas5P(); save_idxs = [2, 1]); verbose = false
)

# Without a symbolic system, a solution saved with `save_idxs` permuting all unknowns is
# detected by the consistency check of the first saved state.
probperm = ODEProblem((u, p, t) -> [-u[1]^3, -2u[2]^3], [1.0, 10.0], (0.0, 0.01))
@test_throws ArgumentError analyze_solution_jacobian(
    solve(probperm, Rodas5P(); save_idxs = [2, 1]); verbose = false
)
@test analyze_solution_jacobian(solve(probperm, Rodas5P()); verbose = false).max_unstable == 0

# Without dense output (`saveat`), the permutation is detected by comparing the first saved
# state with the initial state of the problem.
@test_throws ArgumentError analyze_solution_jacobian(
    solve(probperm, Rodas5P(); save_idxs = [2, 1], saveat = 0.001); verbose = false
)
# With equal initial values, the permutation agrees with the right-hand side at the first
# point and is detected at later saved points.
probperm2 = ODEProblem((u, p, t) -> [-u[1]^3, -(1 + 100t) * u[2]^3], [1.0, 1.0], (0.0, 0.1))
@test_throws ArgumentError analyze_solution_jacobian(
    solve(probperm2, Rodas5P(); save_idxs = [2, 1]); verbose = false
)
@test analyze_solution_jacobian(solve(probperm2, Rodas5P()); verbose = false).max_unstable == 0

# In-place right-hand side under the default `AutoSpecialize` (function wrappers).
probiip = ODEProblem((du, u, p, t) -> (du .= -2 .* u), [1.0, 2.0], (0.0, 1.0))
resiip = analyze_solution_jacobian(solve(probiip, Rodas5P()); verbose = false)
@test all(l -> all(≈(-2), l), resiip.eigenvalues)

# Solutions without dense output (`saveat`) are not checked against the interpolated
# derivative, which is a secant slope there.
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
@test_throws DimensionMismatch analyze_solution_jacobian(solbad; verbose = false, check_state = false)

# The report is produced through the `solution_jacobian_analysis` toggle.
@test_logs (:info, r"Solution Jacobian mode analysis") match_mode = :any analyze_solution_jacobian(sol)
@test_logs min_level = Logging.Info analyze_solution_jacobian(
    sol; verbose = MTKVerbosity(solution_jacobian_analysis = SciMLLogging.Silent)
)
