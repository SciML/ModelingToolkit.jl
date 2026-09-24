using ModelingToolkitBase
using SciMLBase
# narrow import: a blanket `using SymbolicIndexingInterface` would clash with
# ModelingToolkitBase's exported `observed`
using SymbolicIndexingInterface: state_values, parameter_values
using NonlinearSolve   # HomotopySweep (re-exported) + inner nonlinear solvers
using Test

# OMC parity fixture: the Modelica Buildings PressureDrop pattern
#
#     m_flow ~ homotopy(sign(dp)*sqrt(abs(dp))*k, m_flow_nominal*dp/dp_nominal)
#
# Written in `var ~ homotopy(...)` form ON PURPOSE: with the variable itself as
# the LHS, `mtkcompile` eliminates `m_flow` into `observed` (a `0 ~ m_flow - ...`
# form would keep it as an unknown). That makes this the suite's ONLY coverage
# where the homotopy node reaches the swept residual through an OBSERVED equation
# — exercising the observed-rewriting in `lower_homotopy` and the
# observed-inlining in `generate_homotopy_residual` on a real Modelica-parity
# model. The single unknown `w` is pinned THROUGH the homotopy-observed `m_flow`
# by `w*|w| ~ m_flow^2/5`: `u*|u|` is monotone, so the root is unique, and at the
# actual (λ = 1) endpoint `m_flow = √5` gives `w = 1`, the exact OMC solution
# point (`w = 1`, `m_flow = √5`, `dp = 5`).
#
# Frozen reference (MTK CI does not run OMC): OMC v1.27, captured 2026-06-06 from
# ~/code/modelica-OMEdit/homotopy-tests/PressureDropParityFixture.mo —
# m_flow = sqrt(5) = 2.2360679774997896, at dp = 5, k = 1, m_flow_nominal = 1,
# dp_nominal = 5.
@testset "homotopy OMC parity — Buildings PressureDrop fixture (m_flow eliminated into observed)" begin
    OMC_M_FLOW = 2.2360679774997896   # = √5, OMC v1.27 ref captured 2026-06-06
    OMC_W = 1.0
    PARITY_TOL = 1.0e-6

    @variables m_flow w
    @parameters dp k m_flow_nominal dp_nominal
    basic_flow(dpv, kv) = sign(dpv) * sqrt(abs(dpv)) * kv

    @mtkcompile sys = System(
        [
            m_flow ~ homotopy(basic_flow(dp, k), m_flow_nominal * dp / dp_nominal),
            0 ~ w * abs(w) - m_flow^2 / 5,
        ];
        guesses = [w => 0.5, m_flow => 2.0]
    )

    # The key structural fact: m_flow is eliminated into observed, NOT an unknown
    # — the elimination this fixture exists to guard.
    @test !any(v -> isequal(v, m_flow), unknowns(sys))
    @test any(eq -> isequal(eq.lhs, m_flow), observed(sys))
    # ...and the homotopy node reaches the system through OBSERVED, not equations.
    @test ModelingToolkitBase.has_any_homotopy(sys)
    @test !ModelingToolkitBase.has_homotopy_in_equations(equations(sys))
    @test any(eq -> ModelingToolkitBase.has_homotopy(eq.rhs), observed(sys))

    op = [w => 0.5, dp => 5.0, k => 1.0, m_flow_nominal => 1.0, dp_nominal => 5.0]
    prob = HomotopyProblem(sys, op)
    @test prob isa SciMLBase.HomotopyProblem

    # λ-dependence of the wired residual: homotopy sits ONLY in observed here, so
    # differing endpoints come solely from observed-rewriting (`lower_homotopy`)
    # + observed-inlining (`generate_homotopy_residual`) — the codegen path under
    # guard.
    u = copy(state_values(prob))
    p = parameter_values(prob)
    if SciMLBase.isinplace(prob)
        r0 = similar(u)
        r1 = similar(u)
        prob.f(r0, u, p, 0.0)
        prob.f(r1, u, p, 1.0)
        @test !(r0 ≈ r1)
    else
        @test !(prob.f(u, p, 0.0) ≈ prob.f(u, p, 1.0))
    end

    @testset "default sweep matches OMC" begin
        sol = solve(prob; abstol = 1.0e-12, reltol = 1.0e-12)
        @test SciMLBase.successful_retcode(sol)
        @test abs(sol[m_flow] - OMC_M_FLOW) < PARITY_TOL
        @test abs(sol[w] - OMC_W) < PARITY_TOL
    end

    @testset "explicit HomotopySweep matches OMC" begin
        sol = solve(prob, HomotopySweep(); abstol = 1.0e-12, reltol = 1.0e-12)
        @test SciMLBase.successful_retcode(sol)
        @test abs(sol[m_flow] - OMC_M_FLOW) < PARITY_TOL
        @test abs(sol[w] - OMC_W) < PARITY_TOL
    end
end
