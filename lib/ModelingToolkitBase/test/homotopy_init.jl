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
end
