using ModelingToolkit, Test, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D
using NonlinearSolve
using OrdinaryDiffEqBDF
using SciMLBase

# Robertson in index-1 DAE form: the canonical fully implicit test case.
@parameters k1 = 0.04 k2 = 3.0e7 k3 = 1.0e4
@variables y1(t) = 1.0 y2(t) = 0.0 y3(t) = 0.0
@named rob = System(
    [
        D(y1) ~ -k1 * y1 + k3 * y2 * y3,
        D(y2) ~ k1 * y1 - k2 * y2^2 - k3 * y2 * y3,
        0 ~ y1 + y2 + y3 - 1,
    ], t
)
robsys = mtkcompile(rob)
rob_du0 = [D(y1) => -0.04, D(y2) => 0.04]
rob_tspan = (0.0, 1.0e4)

# `w` cannot be solved for symbolically, so the algebraic equation survives `mtkcompile`
# and the compiled system keeps a singular mass matrix: the semi-explicit DAE case.
@variables x(t) = 1.0 w(t) = 0.5
@named semi = System([D(x) ~ -x + w, 0 ~ w^3 + w - x], t)
semisys = mtkcompile(semi)
semi_du0 = [D(x) => -0.5, D(w) => 0.0]

# `nlstep_data` is a stage-system template, so a stage residual has to be evaluated at
# hand-chosen `γ`/`tmp` values rather than at whatever the problem was built with.
function set_stage!(nlp, nd, γ₁, γ₂, c, outer, inner)
    nd.set_γ_c(nlp, (γ₁, γ₂, 1.0, c))
    nd.set_outer_tmp(nlp, outer)
    nd.set_inner_tmp(nlp, inner)
    return nlp
end

@testset "`nlstep_data` is attached to the `DAEFunction`" begin
    prob = DAEProblem(robsys, rob_du0, rob_tspan; nlstep = true)
    nd = prob.f.nlstep_data
    @test nd isa SciMLBase.ODENLStepData

    # the stage unknowns are a subset of the ODE unknowns: tearing introduced nothing new
    @test all(!isnothing, nd.u0perm)
    @test allunique(nd.u0perm)
    @test issubset(nd.u0perm, eachindex(unknowns(robsys)))
    @test length(nd.nlprob.u0) == length(nd.u0perm)

    # γ₃ is unused by the fully implicit stage system but must stay settable, so that
    # `set_γ_c` takes the same four values as it does on the mass-matrix path
    nlp = deepcopy(nd.nlprob)
    nd.set_γ_c(nlp, (2.0, 3.0, 1.0, 0.5))

    @test DAEProblem(robsys, rob_du0, rob_tspan).f.nlstep_data === nothing
end

# `g(z, p') = F(γ₁ z + outer_tmp, γ₂ z + inner_tmp, p, c)`. With `nlstep_compile = false`
# the stage system is not teared, so its residual is comparable to `F` componentwise.
@testset "stage residual is the fully implicit residual: $name" for (name, sys, du0) in (
        ("Robertson", robsys, rob_du0), ("semi-explicit", semisys, semi_du0),
    )
    prob = DAEProblem(sys, du0, (0.0, 1.0); nlstep = true, nlstep_compile = false)
    nd = prob.f.nlstep_data
    N = length(unknowns(sys))
    @test nd.u0perm == 1:N

    γ₁, γ₂, c = 7.0, 3.0, 0.3
    outer = [0.11, -0.22]
    inner = [0.05, 0.07]
    nlp = set_stage!(deepcopy(nd.nlprob), nd, γ₁, γ₂, c, outer, inner)

    for z in ([0.9, 0.01], [0.4, -0.3])
        resid = zeros(N)
        nlp.f(resid, z, nlp.p)
        expected = zeros(N)
        prob.f(expected, γ₁ .* z .+ outer, γ₂ .* z .+ inner, prob.p, c)
        @test resid ≈ expected
    end
end

# The teared stage system drops `x`, so this also covers the case where `u0perm` is a
# strict subset of the ODE unknowns and `nlprobmap` has to rebuild the rest.
@testset "compiled stage solution solves the fully implicit residual: $name" for (
        name, sys, du0, dt, n_stage,
    ) in (
        ("Robertson", robsys, rob_du0, 1.0e-4, 2),
        ("semi-explicit", semisys, semi_du0, 1.0e-2, 1),
    )
    prob = DAEProblem(sys, du0, (0.0, 1.0); nlstep = true)
    nd = prob.f.nlstep_data
    N = length(unknowns(sys))
    # Which unknowns tearing retains is the contract; the order it returns them in is not,
    # and does vary by Julia version (Robertson gives `[1, 2]` on 1.11/1.12 and `[2, 1]` on
    # 1.10). `nlprobmap` undoes the permutation either way, which the residual check below
    # verifies, so pin the count rather than the arrangement.
    @test length(nd.u0perm) == n_stage
    @test allunique(nd.u0perm)
    @test issubset(nd.u0perm, eachindex(unknowns(sys)))

    # a backward Euler step of size `dt` from `uprev`: `du ≈ (u - uprev) / dt`, i.e.
    # `γ₁ = inv(dt)`, `outer_tmp = -uprev / dt`, `γ₂ = 1`, `inner_tmp = 0`
    uprev = prob.u0
    nlp = set_stage!(
        deepcopy(nd.nlprob), nd, inv(dt), 1.0, dt, -uprev ./ dt, zeros(N)
    )

    # `γ₁ z + outer_tmp` is `(z - uprev) / dt`, so the stage residual is a cancellation of
    # two `O(1/dt)` quantities and cannot be resolved below an absolute `eps() / dt`
    # (`2.2e-12` at `dt = 1e-4`). A tolerance under that floor makes Newton report
    # `Stalled` at a fully converged iterate, so keep it above.
    sol = solve(nlp, NewtonRaphson(); abstol = 1.0e-10, reltol = 1.0e-10)
    @test SciMLBase.successful_retcode(sol)

    z = nd.nlprobmap(sol)
    @test length(z) == N
    resid = zeros(N)
    prob.f(resid, (z .- uprev) ./ dt, z, prob.p, dt)
    @test maximum(abs, resid) < 1.0e-8
end

# `γ₁` is the `gamma` of the `DAEFunction` Jacobian signature `jac(J, du, u, p, gamma, t)`:
# the stage residual's Jacobian in `z` is `γ₁ * dF/d(du) + dF/du`.
@testset "`γ₁` matches the `DAEFunction` Jacobian's `gamma`" begin
    prob = DAEProblem(
        semisys, semi_du0, (0.0, 1.0); jac = true, nlstep = true, nlstep_compile = false
    )
    nd = prob.f.nlstep_data
    N = length(unknowns(semisys))
    @test nd.nlprob.f.jac !== nothing

    γ₁, γ₂, c = 5.0, 1.0, 0.25
    outer = [0.3, 0.0]
    nlp = set_stage!(deepcopy(nd.nlprob), nd, γ₁, γ₂, c, outer, zeros(N))

    z = [0.8, 0.6]
    Jstage = zeros(N, N)
    nlp.f.jac(Jstage, z, nlp.p)
    Jdae = zeros(N, N)
    prob.f.jac(Jdae, γ₁ .* z .+ outer, γ₂ .* z, prob.p, γ₁, c)
    @test Jstage ≈ Jdae
end

@testset "attaching `nlstep_data` does not change the solve" begin
    ref = solve(
        DAEProblem(robsys, rob_du0, rob_tspan), DFBDF(); abstol = 1.0e-10, reltol = 1.0e-10
    )
    withnl = solve(
        DAEProblem(robsys, rob_du0, rob_tspan; nlstep = true), DFBDF();
        abstol = 1.0e-10, reltol = 1.0e-10
    )
    @test SciMLBase.successful_retcode(withnl)
    @test withnl.t ≈ ref.t
    @test withnl.u ≈ ref.u
end

@testset "`nlstep_scc` builds an SCC stage problem" begin
    prob = DAEProblem(robsys, rob_du0, rob_tspan; nlstep = true, nlstep_scc = true)
    nd = prob.f.nlstep_data
    @test nd.nlprob isa SciMLBase.AbstractNonlinearProblem
    @test all(!isnothing, nd.u0perm)
end
