using Test
using ModelingToolkitBase
using ModelingToolkitBase: System, mtkcompile, HomotopySweep, TrivialHomotopy
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using OrdinaryDiffEqRosenbrock: Rodas5P
using OrdinaryDiffEqNonlinearSolve
using SciMLBase
using NonlinearSolve: NewtonRaphson

@testset "homotopy OMC parity — Buildings PressureDrop fixture" begin
    # Reference values captured 2026-05-26 from OMC v1.27.0-dev running
    # `~/code/modelica-OMEdit/homotopy-tests/PressureDropParityFixture.mo`
    # via the sandbox's `run_mo.sh`. m_flow(0) is √5 (the actual root);
    # x(0) is the user-supplied start value, untouched by init.
    OMC_M_FLOW_AT_T0 = 2.2360679774997896   # = sqrt(5)
    OMC_X_AT_T0      = 1.0
    PARITY_TOL       = 1e-6

    @variables x(t) m_flow(t)
    @parameters dp k m_flow_nominal dp_nominal x_decay

    basic_flow(dp_val, k_val) = sign(dp_val) * sqrt(abs(dp_val)) * k_val

    eqs = [D(x) ~ -x_decay * x,
           0 ~ m_flow - homotopy(
                   basic_flow(dp, k),
                   m_flow_nominal * dp / dp_nominal)]
    @named sys = System(eqs, t; guesses = [m_flow => 2.0])
    sys = mtkcompile(sys)
    op = Dict(x => 1.0, dp => 5.0, k => 1.0,
              m_flow_nominal => 1.0, dp_nominal => 5.0, x_decay => 1.0)
    prob = ODEProblem(sys, op, (0.0, 1.0))

    @testset "default (TrivialThenSweep) matches OMC" begin
        sol = solve(prob, Rodas5P(); abstol = 1e-12, reltol = 1e-12)
        @test SciMLBase.successful_retcode(sol)
        @test abs(sol[m_flow][1] - OMC_M_FLOW_AT_T0) < PARITY_TOL
        @test abs(sol[x][1]      - OMC_X_AT_T0)      < PARITY_TOL
    end

    @testset "explicit HomotopySweep matches OMC" begin
        meta = prob.f.initialization_data.metadata
        sweep = HomotopySweep(; inner = NewtonRaphson(),
                                schedule = 0.0:0.1:1.0,
                                set_λ! = meta.homotopy_set_λ!)
        initalg = SciMLBase.OverrideInit(; nlsolve = sweep)
        sol = solve(prob, Rodas5P(); initializealg = initalg,
                    abstol = 1e-12, reltol = 1e-12)
        @test SciMLBase.successful_retcode(sol)
        @test abs(sol[m_flow][1] - OMC_M_FLOW_AT_T0) < PARITY_TOL
    end

    @testset "explicit TrivialHomotopy matches OMC" begin
        initalg = SciMLBase.OverrideInit(; nlsolve = TrivialHomotopy())
        sol = solve(prob, Rodas5P(); initializealg = initalg,
                    abstol = 1e-12, reltol = 1e-12)
        @test SciMLBase.successful_retcode(sol)
        @test abs(sol[m_flow][1] - OMC_M_FLOW_AT_T0) < PARITY_TOL
    end
end
