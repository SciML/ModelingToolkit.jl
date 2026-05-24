using Test
using ModelingToolkitBase
using ModelingToolkitBase: rewrite_with_lambda, rewrite_with_lambda_in_equations,
                           has_homotopy, has_homotopy_in_equations
using Symbolics

@testset "homotopy operator — L1 lowering rewrite_with_lambda" begin
    @testset "L1-Q1 basic single node" begin
        @variables x
        @parameters p
        expr = homotopy(x^2 - p, x - sqrt(p))
        rewritten, λ = rewrite_with_lambda(expr)
        @test !has_homotopy(rewritten)
        at_one  = Symbolics.simplify(Symbolics.substitute(rewritten, Dict(λ => 1.0)))
        at_zero = Symbolics.simplify(Symbolics.substitute(rewritten, Dict(λ => 0.0)))
        # At λ=1 the lowered expression reduces to `actual`; at λ=0 to `simplified`.
        @test isequal(Symbolics.unwrap(at_one),  Symbolics.unwrap(x^2 - p))
        @test isequal(Symbolics.unwrap(at_zero), Symbolics.unwrap(x - sqrt(p)))
    end

    @testset "L1-Q2 nested homotopy collapses with single λ" begin
        @variables x
        @parameters p
        inner = homotopy(x^2 - p, x - sqrt(p))
        outer = homotopy(inner, x - 1)
        rewritten, λ = rewrite_with_lambda(outer)
        @test !has_homotopy(rewritten)
        at_one = Symbolics.simplify(Symbolics.substitute(rewritten, Dict(λ => 1.0)))
        # Nested homotopy at λ=1 collapses fully to the innermost `actual`.
        @test isequal(Symbolics.unwrap(at_one), Symbolics.unwrap(x^2 - p))
    end

    @testset "L1-Q3 multiple nodes share the same λ" begin
        @variables x y
        @parameters p q
        e1 = homotopy(x^2 - p, x - 1)
        e2 = homotopy(y^2 - q, y - 1)
        eqs = [0 ~ e1, 0 ~ e2]
        new_eqs, λ = rewrite_with_lambda_in_equations(eqs)
        @test length(new_eqs) == 2
        params_in_rhs = Set{Symbol}()
        for eq in new_eqs
            for var in Symbolics.get_variables(eq.rhs)
                push!(params_in_rhs, nameof(Symbolics.unwrap(var)))
            end
        end
        @test :__homotopy_λ in params_in_rhs
        @test !has_homotopy_in_equations(new_eqs)
    end

    @testset "L1-Q4 no-homotopy equations pass through unchanged" begin
        @variables x
        @parameters p
        eqs = [0 ~ x^2 - p, 0 ~ x + 1]
        new_eqs, λ = rewrite_with_lambda_in_equations(eqs)
        @test isequal(new_eqs[1].rhs, eqs[1].rhs)
        @test isequal(new_eqs[2].rhs, eqs[2].rhs)
        @test λ !== nothing
    end
end
