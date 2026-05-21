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
end
