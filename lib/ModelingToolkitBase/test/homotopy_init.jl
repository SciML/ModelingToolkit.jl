using Test
using ModelingToolkitBase
using ModelingToolkitBase: rewrite_trivial, has_homotopy
using Symbolics

@testset "homotopy operator — L0 trivial rewrite" begin
    @testset "Q1 basic" begin
        @variables x
        @parameters p
        expr = homotopy(x^2 - p, x - sqrt(p))
        rewritten = rewrite_trivial(expr)
        @test isequal(Symbolics.unwrap(rewritten), Symbolics.unwrap(x^2 - p))
        @test !has_homotopy(rewritten)
    end

    @testset "Q2 nested" begin
        @variables x
        @parameters p
        inner = homotopy(x^2 - p, x - sqrt(p))
        outer = homotopy(inner, x - 1)
        @test isequal(Symbolics.unwrap(rewrite_trivial(outer)),
                      Symbolics.unwrap(x^2 - p))

        triple = homotopy(homotopy(homotopy(x, x + 1), x + 2), x + 3)
        @test isequal(Symbolics.unwrap(rewrite_trivial(triple)),
                      Symbolics.unwrap(x))

        @test !has_homotopy(rewrite_trivial(outer))
        @test !has_homotopy(rewrite_trivial(triple))
    end

    @testset "Q3 Base.ifelse branches" begin
        @variables x
        @parameters p
        cond = p > 0
        branch_expr = Base.ifelse(
            cond,
            homotopy(x^2 - p, x - sqrt(abs(p))),
            homotopy(-(x^2) - p, x + sqrt(abs(p))),
        )
        rewritten = rewrite_trivial(branch_expr)

        # Outer ifelse structure preserved
        @test occursin("ifelse", repr(rewritten))
        # Both branches' homotopy nodes are gone
        @test !has_homotopy(rewritten)

        # Folding the cond to true/false picks the right actual
        true_folded  = Symbolics.simplify(Symbolics.substitute(rewritten, Dict(cond => true)))
        false_folded = Symbolics.simplify(Symbolics.substitute(rewritten, Dict(cond => false)))
        @test isequal(Symbolics.unwrap(true_folded),  Symbolics.unwrap(x^2 - p))
        @test isequal(Symbolics.unwrap(false_folded), Symbolics.unwrap(-(x^2) - p))
    end

    @testset "Q4 broadcast" begin
        @variables x[1:3] p[1:3]
        actuals     = [x[i]^2 - p[i] for i in 1:3]
        simplifieds = [x[i]   - p[i] for i in 1:3]
        vec_expr  = homotopy.(actuals, simplifieds)
        rewritten = rewrite_trivial.(vec_expr)
        @test length(rewritten) == 3
        for i in 1:3
            @test isequal(Symbolics.unwrap(rewritten[i]),
                          Symbolics.unwrap(actuals[i]))
            @test !has_homotopy(rewritten[i])
        end
    end

    @testset "Integration: homotopy in NonlinearSystem init" begin
        using ModelingToolkitBase: System, mtkcompile
        using ModelingToolkitBase: InitializationProblem
        using NonlinearSolve: NewtonRaphson, solve
        using SciMLBase

        @variables x
        @parameters p
        # Equation: homotopy(actual = x^2 - p, simplified = x - 1) = 0
        # L0 trivial must solve actual: x^2 = p ⇒ x = ±√p.
        # If rewrite is broken and simplified got selected: x = 1.
        # At p = 9: actual root x = 3 (or -3); simplified root x = 1.
        # The x^2 ≈ p assertion distinguishes the two outcomes.
        eqs = [0 ~ homotopy(x^2 - p, x - 1)]
        @named sys = System(eqs, [x], [p])
        sys = mtkcompile(sys)

        prob = InitializationProblem{false}(sys, nothing, Dict(p => 9.0, x => 2.5))
        sol = solve(prob, NewtonRaphson(); abstol = 1e-10, reltol = 1e-10)
        @test SciMLBase.successful_retcode(sol)
        @test abs(sol.u[1]^2 - 9.0) < 1e-6   # actual equation solved (NOT simplified, which would give x=1)
        @test abs(sol.u[1] - 1.0) > 0.5      # not the simplified root
    end
end
