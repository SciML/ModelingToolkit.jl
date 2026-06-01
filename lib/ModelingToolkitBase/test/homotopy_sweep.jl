using Test
using ModelingToolkitBase
using ModelingToolkitBase: HomotopySweep, TrivialHomotopy, TrivialThenSweep
using SciMLBase
using NonlinearSolve: NewtonRaphson

@testset "Homotopy nlsolve algorithms" begin
    # Build the canonical homotopy fixture: f(u, p) = (1-λ)*(u-1) + λ*(u^2-4)
    # - At λ=0: u = 1 (simplified).  At λ=1: u = ±2 (actual).
    # - From u0 = 1, sweep should track to u ≈ 2.
    f_canonical! = (du, u, p) -> begin
        λ = p[1]
        du[1] = (1 - λ) * (u[1] - 1) + λ * (u[1]^2 - 4)
    end
    set_λ_explicit = (p, v) -> begin
        p2 = copy(p)
        p2[1] = v
        p2
    end

    @testset "S1 HomotopySweep linear schedule converges to actual" begin
        prob = NonlinearProblem(f_canonical!, [1.0], [0.0])
        alg = HomotopySweep(; inner = NewtonRaphson(),
                              schedule = 0.0:0.1:1.0,
                              set_λ! = set_λ_explicit)
        sol = solve(prob, alg)
        @test SciMLBase.successful_retcode(sol)
        @test abs(sol.u[1] - 2.0) < 1e-6
    end

    @testset "S2 HomotopySweep step failure returns unsuccessful retcode" begin
        # Force impossible homotopy: f = u^2 + 1, no real root → Newton diverges.
        f_bad! = (du, u, p) -> begin
            λ = p[1]
            du[1] = (1 - λ) * u[1] + λ * (u[1]^2 + 1)
        end
        prob = NonlinearProblem(f_bad!, [0.5], [0.0])
        alg = HomotopySweep(; inner = NewtonRaphson(),
                              schedule = 0.0:0.1:1.0,
                              set_λ! = set_λ_explicit,
                              maxiters_per_step = 20)
        sol = solve(prob, alg)
        @test !SciMLBase.successful_retcode(sol)
    end

    @testset "S2b failed step resets λ to actual form (1.0) before returning" begin
        # SWEEP-2: when a sweep step fails mid-schedule, the returned solution's
        # parameter vector must hold λ = 1.0 (actual form), not the intermediate
        # λ where continuation stalled. Otherwise a downstream consumer reading
        # the failed solution's parameters sees a half-homotopied system.
        # Fixture: f = (1-λ)*u + λ*(u^2+1). A real root exists only for small λ;
        # around λ ≈ 0.5 the quadratic loses its real root and Newton diverges,
        # so the failing step's λ is strictly between 0 and 1.
        f_bad! = (du, u, p) -> begin
            λ = p[1]
            du[1] = (1 - λ) * u[1] + λ * (u[1]^2 + 1)
        end
        prob = NonlinearProblem(f_bad!, [0.5], [0.0])
        alg = HomotopySweep(; inner = NewtonRaphson(),
                              schedule = 0.0:0.1:1.0,
                              set_λ! = set_λ_explicit,
                              maxiters_per_step = 20)
        sol = solve(prob, alg)
        @test !SciMLBase.successful_retcode(sol)
        @test sol.prob.p[1] == 1.0
    end

    @testset "S3 HomotopySweep default constructor" begin
        alg = HomotopySweep(; set_λ! = ((p, v) -> p))
        @test alg.schedule == 0.0:0.1:1.0
        # NonlinearSolve 5+ exposes `NewtonRaphson` as a constructor function
        # (not a type), so compare against the type of an instance.
        @test alg.inner isa typeof(NewtonRaphson())
    end

    @testset "S3b default inner supplied by NonlinearSolve extension" begin
        # MERGE-3: the default inner solver is provided by the
        # `MTKNonlinearSolveExt` package extension (populated when NonlinearSolve
        # is loaded), not by an eager construction-time `Base.require`. With
        # NonlinearSolve loaded the factory Ref is set and resolves to
        # `NonlinearSolve.NewtonRaphson`; without it, `_default_inner` falls back
        # to `SimpleNonlinearSolve.SimpleNewtonRaphson` (a hard dependency).
        @test ModelingToolkitBase._DEFAULT_INNER_FACTORY[] !== nothing
        @test ModelingToolkitBase._default_inner() isa typeof(NewtonRaphson())
    end

    @testset "S4 TrivialHomotopy dispatches to inner once" begin
        f! = (du, u, p) -> du[1] = u[1]^2 - 4
        prob = NonlinearProblem(f!, [1.5], Float64[])
        alg = TrivialHomotopy()
        sol = solve(prob, alg)
        @test SciMLBase.successful_retcode(sol)
        @test abs(sol.u[1] - 2.0) < 1e-6
    end

    @testset "S5 TrivialThenSweep succeeds via trivial path when guess is good" begin
        # u0 = 1.9 is close to actual root u=2 → Newton converges from trivial.
        # Sweep should NOT be invoked. We detect by `sol.original.path === :trivial`.
        # NOTE: p[1] = 1.0 is intentional — TrivialHomotopy does NOT touch λ,
        # so the canonical fixture must already be in actual-form. Don't "fix"
        # this to 0.0; that would silently flip the trivial path to simplified.
        prob = NonlinearProblem(f_canonical!, [1.9], [1.0])  # p[1]=1.0 → already trivial
        alg = TrivialThenSweep(;
            trivial = TrivialHomotopy(),
            sweep = HomotopySweep(; inner = NewtonRaphson(),
                                    schedule = 0.0:0.1:1.0,
                                    set_λ! = set_λ_explicit))
        sol = solve(prob, alg)
        @test SciMLBase.successful_retcode(sol)
        @test abs(sol.u[1] - 2.0) < 1e-6
        @test sol.original isa NamedTuple && sol.original.path === :trivial
    end

    @testset "S6 TrivialThenSweep falls back to sweep when trivial fails" begin
        # Fixture: simplified u (linear, root u=0), actual atan(u) (root u=0).
        # Newton on atan(u) from u0=10 is the classic diverging case (stalls/diverges
        # because |u_{n+1}| > |u_n| outside the basin of attraction).
        # Sweep from u0=10 at λ=0 immediately drives u→0 via the simplified branch,
        # then atan walks smoothly to root u=0 across the schedule.
        f_hard! = (du, u, p) -> begin
            λ = p[1]
            du[1] = (1 - λ) * u[1] + λ * atan(u[1])
        end
        prob = NonlinearProblem(f_hard!, [10.0], [1.0])  # start at λ=1 (trivial)
        alg = TrivialThenSweep(;
            trivial = TrivialHomotopy(; inner = NewtonRaphson()),
            sweep = HomotopySweep(; inner = NewtonRaphson(),
                                    schedule = 0.0:0.1:1.0,
                                    set_λ! = set_λ_explicit))
        sol = solve(prob, alg)
        @test SciMLBase.successful_retcode(sol)
        @test abs(atan(sol.u[1])) < 1e-4
        @test sol.original.path === :sweep_fallback
    end
end
