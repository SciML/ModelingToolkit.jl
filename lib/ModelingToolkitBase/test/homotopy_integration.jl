using Test
using ModelingToolkitBase
using ModelingToolkitBase: System, mtkcompile,
                           HomotopySweep, TrivialHomotopy, TrivialThenSweep
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using OrdinaryDiffEqRosenbrock: Rodas5P
using OrdinaryDiffEqNonlinearSolve
using SciMLBase
using NonlinearSolve: NewtonRaphson

# Shared fixture: algebraic 0 = homotopy(y^2 - p, y - 1).
# Simplified at λ=0: y = 1.  Actual at λ=1: y = ±√p.  At p=9, guess y=2.5 → y=3.
function build_homotopy_fixture(; guess = 2.5, p_val = 9.0)
    @variables x(t) y(t)
    @parameters p
    eqs = [D(x) ~ -x,
           0 ~ homotopy(y^2 - p, y - 1)]
    @named sys = System(eqs, t; guesses = [y => guess])
    sys = mtkcompile(sys)
    prob = ODEProblem(sys, Dict(x => 1.0, p => p_val), (0.0, 1.0))
    return prob, y
end

@testset "homotopy operator — integration end-to-end" begin
    @testset "I1 default (TrivialThenSweep) on good guess → trivial path wins" begin
        prob, y = build_homotopy_fixture(; guess = 2.5)
        sol = solve(prob, Rodas5P(); abstol = 1e-10, reltol = 1e-10)
        @test SciMLBase.successful_retcode(sol)
        @test abs(sol[y][1] - 3.0) < 1e-4
    end

    @testset "I2 explicit HomotopySweep" begin
        prob, y = build_homotopy_fixture(; guess = 2.5)
        meta = prob.f.initialization_data.metadata
        sweep = HomotopySweep(; inner = NewtonRaphson(),
                                schedule = 0.0:0.1:1.0,
                                set_λ! = meta.homotopy_set_λ!)
        initalg = SciMLBase.OverrideInit(; nlsolve = sweep)
        sol = solve(prob, Rodas5P(); initializealg = initalg,
                    abstol = 1e-10, reltol = 1e-10)
        @test SciMLBase.successful_retcode(sol)
        @test abs(sol[y][1] - 3.0) < 1e-4
    end

    @testset "I3 explicit TrivialHomotopy" begin
        prob, y = build_homotopy_fixture(; guess = 2.5)
        initalg = SciMLBase.OverrideInit(; nlsolve = TrivialHomotopy())
        sol = solve(prob, Rodas5P(); initializealg = initalg,
                    abstol = 1e-10, reltol = 1e-10)
        @test SciMLBase.successful_retcode(sol)
        @test abs(sol[y][1] - 3.0) < 1e-4
    end

    @testset "I4 default-injection wires TrivialThenSweep all the way to solve" begin
        # Structural assertion: the OMC-aligned default `TrivialThenSweep`
        # really reaches `prob.kwargs.initializealg` so `solve` consumes it
        # without explicit user opt-in. Trivial-vs-sweep dispatch within
        # `TrivialThenSweep.solve` is tested directly in `homotopy_sweep.jl`
        # (S5 trivial-success, S6 sweep-fallback) — no need to re-prove it
        # through the ODEProblem pipeline.
        prob, y = build_homotopy_fixture(; guess = 2.5)
        @test haskey(prob.kwargs, :initializealg)
        injected = prob.kwargs[:initializealg]
        @test injected isa SciMLBase.OverrideInit
        @test injected.nlsolve isa TrivialThenSweep
        @test injected.nlsolve.trivial isa TrivialHomotopy
        @test injected.nlsolve.sweep isa HomotopySweep
        # And the actual end-to-end solve must succeed using only this default
        sol = solve(prob, Rodas5P(); abstol = 1e-10, reltol = 1e-10)
        @test SciMLBase.successful_retcode(sol)
        @test abs(sol[y][1] - 3.0) < 1e-4
    end
end
