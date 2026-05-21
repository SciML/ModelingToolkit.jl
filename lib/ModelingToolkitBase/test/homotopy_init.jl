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
end
