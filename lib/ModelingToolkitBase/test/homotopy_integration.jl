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

    @testset "I5 out-of-basin guess: continuation rescues init the trivial solve can't" begin
        # The sweep's reason-for-being, exercised THROUGH the real ODEProblem
        # pipeline (homotopy_sweep.jl S6 only proves it for a hand-built
        # NonlinearProblem). `actual = atan(y - 3)` has its single root at y = 3
        # but Newton diverges from any guess outside |y-3| ≲ 1.39 (the classic
        # atan basin-escape). `simplified = y` has root y = 0. The two roots
        # differ, so the landing value distinguishes "continuation tracked
        # λ:0→1 onto the actual root" (y = 3) from "stuck at simplified" (y = 0).
        @variables x(t) y(t)
        eqs = [D(x) ~ -x,
               0 ~ homotopy(atan(y - 3.0), y)]
        # guess y = 12 is far outside the actual root's Newton basin: a single
        # Newton at λ=1 diverges, so reaching y = 3 is only possible via the
        # default TrivialThenSweep continuation walking λ from 0 to 1.
        @named sys = System(eqs, t; guesses = [y => 12.0])
        sys = mtkcompile(sys)

        prob = ODEProblem(sys, Dict(x => 1.0), (0.0, 1.0))
        sol = solve(prob, Rodas5P(); abstol = 1e-10, reltol = 1e-10)
        @test SciMLBase.successful_retcode(sol)
        @test abs(sol[y][1] - 3.0) < 1e-6   # actual root via the sweep
        @test abs(sol[y][1]) > 0.5          # definitively NOT the simplified root y = 0
    end

    @testset "I6 explicit TrivialHomotopy on an unsolvable init fails cleanly" begin
        # Regression: an init the single-Newton trivial path cannot solve must
        # surface as an unsuccessful retcode, NOT a cryptic internal error. The
        # DAE-init driver forwards its `initializealg` (an OverrideInit) into the
        # nlsolve `solve`; if that leaks into the inner solver it re-runs the
        # OverrideInit on a problem with no initialization_data and throws. Same
        # out-of-basin fixture as I5, but pinning the trivial path (no sweep
        # fallback) so init genuinely cannot converge.
        @variables x(t) y(t)
        eqs = [D(x) ~ -x,
               0 ~ homotopy(atan(y - 3.0), y)]
        @named sys = System(eqs, t; guesses = [y => 12.0])
        sys = mtkcompile(sys)

        prob = ODEProblem(sys, Dict(x => 1.0), (0.0, 1.0);
                          initializealg = SciMLBase.OverrideInit(nlsolve = TrivialHomotopy()))
        sol = solve(prob, Rodas5P())                  # must NOT throw an internal error
        @test !SciMLBase.successful_retcode(sol)      # clean init failure instead
    end
end
