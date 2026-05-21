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

    @testset "Integration: homotopy in ODESystem init" begin
        using ModelingToolkitBase: System, mtkcompile
        using ModelingToolkitBase: InitializationProblem
        using ModelingToolkitBase: t_nounits as t, D_nounits as D
        using ModelingToolkitBase: has_homotopy_in_equations, equations
        using NonlinearSolve: NewtonRaphson, solve
        using SciMLBase

        @variables x(t) y(t)
        @parameters p
        # `y` is algebraic, constrained by homotopy(actual = y^2 - p, simplified = y - 1).
        # L0 trivial rewrite must replace the constraint with `y^2 - p = 0`, so init solves
        # to y = ±√p. At p = 9 with guess y = 2.5, init converges to y ≈ 3. If the rewrite
        # never fired, the init system would carry the opaque `homotopy(...)` symbolic call
        # — which is the structural assertion below.
        eqs = [D(x) ~ -x,
               0  ~ homotopy(y^2 - p, y - 1)]
        @named sys = System(eqs, t; guesses = [y => 2.5])
        sys = mtkcompile(sys)

        prob = InitializationProblem{false}(sys, 0.0, Dict(x => 1.0, p => 9.0))

        # Structural: after the init pipeline applies rewrite_trivial, no equation in
        # the wrapped initialization system should still contain a `homotopy(...)` node.
        @test !has_homotopy_in_equations(equations(prob.f.sys))

        # Numerical: init must converge to the actual root.
        sol = solve(prob, NewtonRaphson(); abstol = 1e-10, reltol = 1e-10)
        @test SciMLBase.successful_retcode(sol)
        @test abs(sol.u[1]^2 - 9.0) < 1e-6   # actual equation y^2 = p satisfied
        @test abs(sol.u[1] - 1.0) > 0.5      # not the simplified root (y = 1)
    end
end
