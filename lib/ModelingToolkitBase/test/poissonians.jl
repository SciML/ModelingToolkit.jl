using ModelingToolkitBase, Test
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using ModelingToolkitBase: ispoissonian, getpoissonianrate, POISSONIAN, getvariabletype
using Symbolics: unwrap

# ============================================================================
# Phase 1: Macro and Metadata Tests
# ============================================================================

@testset "Poissonian Macro and Metadata" begin
    @testset "Single poissonian declaration" begin
        @parameters λ
        @poissonians dN(λ)

        # Check that it's marked as a poissonian
        @test ispoissonian(dN)
        @test getvariabletype(unwrap(dN)) === POISSONIAN

        # Check that the rate is stored correctly
        @test isequal(getpoissonianrate(dN), λ)
    end

    @testset "Multiple inline poissonians" begin
        @parameters λ₁ λ₂
        @poissonians dN₁(λ₁) dN₂(λ₂)

        @test ispoissonian(dN₁)
        @test ispoissonian(dN₂)
        @test isequal(getpoissonianrate(dN₁), λ₁)
        @test isequal(getpoissonianrate(dN₂), λ₂)
    end

    @testset "Block syntax poissonians" begin
        @parameters β γ
        @variables S(t) I(t)

        @poissonians begin
            dN_inf(β * S * I)
            dN_rec(γ * I)
        end

        @test ispoissonian(dN_inf)
        @test ispoissonian(dN_rec)
        @test isequal(getpoissonianrate(dN_inf), β * S * I)
        @test isequal(getpoissonianrate(dN_rec), γ * I)
    end

    @testset "State-dependent rates" begin
        @parameters β γ
        @variables S(t) I(t) R(t)

        @poissonians dN(β * S * I)

        # Rate should be the full expression
        rate = getpoissonianrate(dN)
        @test isequal(rate, β * S * I)
    end

    @testset "Constant rate" begin
        @parameters λ
        @poissonians dN(λ)

        @test isequal(getpoissonianrate(dN), λ)
    end

    @testset "Non-poissonian returns nothing for rate" begin
        @variables x(t)

        # Regular variables should return nothing for rate
        @test getpoissonianrate(x) === nothing
        @test !ispoissonian(x)
    end

    @testset "Error on missing rate" begin
        # This should error - rate is required
        @test_throws Exception @macroexpand @poissonians dN
    end
end

# ============================================================================
# Phase 2: System Storage and Accessors Tests
# ============================================================================

@testset "System Storage of Poissonians" begin
    @testset "Explicit constructor with poissonians kwarg" begin
        @parameters β γ
        @variables S(t) I(t) R(t)
        @poissonians dN_inf(β * S * I) dN_rec(γ * I)

        eqs = [
            D(S) ~ -dN_inf,
            D(I) ~ dN_inf - dN_rec,
            D(R) ~ dN_rec
        ]

        # Using explicit 5-argument constructor with poissonians kwarg
        @named sys = System(eqs, t, [S, I, R], [β, γ]; poissonians = [dN_inf, dN_rec])

        @test length(ModelingToolkitBase.poissonians(sys)) == 2
        @test Set(ModelingToolkitBase.poissonians(sys)) == Set([dN_inf, dN_rec])
    end

    @testset "Two-argument constructor auto-detects poissonians" begin
        @parameters β γ
        @variables S(t) I(t) R(t)
        @poissonians dN_inf(β * S * I) dN_rec(γ * I)

        eqs = [
            D(S) ~ -dN_inf,
            D(I) ~ dN_inf - dN_rec,
            D(R) ~ dN_rec
        ]

        # Two-argument constructor should auto-detect poissonians
        @named sys = System(eqs, t)

        @test length(ModelingToolkitBase.poissonians(sys)) == 2
        @test Set(ModelingToolkitBase.poissonians(sys)) == Set([dN_inf, dN_rec])
    end

    @testset "Variables from rate expressions are extracted" begin
        @parameters β γ
        @variables S(t) I(t) R(t)
        @poissonians dN_inf(β * S * I) dN_rec(γ * I)

        eqs = [
            D(S) ~ -dN_inf,
            D(I) ~ dN_inf - dN_rec,
            D(R) ~ dN_rec
        ]

        @named sys = System(eqs, t)

        # S, I, R should be in unknowns (extracted from rate expressions and equations)
        sys_unknowns = Set(ModelingToolkitBase.unknowns(sys))
        @test S in sys_unknowns
        @test I in sys_unknowns
        @test R in sys_unknowns

        # β, γ should be in parameters (extracted from rate expressions)
        sys_params = Set(ModelingToolkitBase.parameters(sys))
        @test β in sys_params
        @test γ in sys_params
    end

    @testset "Poissonians are not included in unknowns" begin
        @parameters λ
        @variables X(t)
        @poissonians dN(λ)

        eqs = [D(X) ~ dN]

        @named sys = System(eqs, t)

        # dN should not be in unknowns
        @test !any(x -> isequal(x, dN), ModelingToolkitBase.unknowns(sys))
        # dN should be in poissonians
        @test any(x -> isequal(x, dN), ModelingToolkitBase.poissonians(sys))
    end

    @testset "get_poissonians accessor" begin
        @parameters λ
        @variables X(t)
        @poissonians dN(λ)

        eqs = [D(X) ~ dN]
        @named sys = System(eqs, t)

        # Test the low-level accessor
        @test length(ModelingToolkitBase.get_poissonians(sys)) == 1
        @test any(x -> isequal(x, dN), ModelingToolkitBase.get_poissonians(sys))
    end

    @testset "Empty poissonians for systems without them" begin
        @parameters a
        @variables X(t)

        eqs = [D(X) ~ -a * X]
        @named sys = System(eqs, t)

        @test isempty(ModelingToolkitBase.poissonians(sys))
    end
end
