using ModelingToolkitBase, Test, JumpProcesses
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using ModelingToolkitBase: ispoissonian, getpoissonianrate, POISSONIAN, getvariabletype
using JumpProcesses: ConstantRateJump, VariableRateJump
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

# ============================================================================
# Phase 3: Poissonian to Jump Conversion Tests
# ============================================================================

@testset "Poissonian to Jump Conversion" begin
    @testset "Pure-jump SIR model - VariableRateJumps" begin
        @parameters β γ
        @variables S(t) I(t) R(t)
        @poissonians dN_inf(β * S * I) dN_rec(γ * I)

        eqs = [
            D(S) ~ -dN_inf,
            D(I) ~ dN_inf - dN_rec,
            D(R) ~ dN_rec
        ]

        @named sys = System(eqs, t)
        compiled_sys = mtkcompile(sys)

        # All equations should be removed (pure-jump system)
        @test isempty(equations(compiled_sys))

        # Should have 2 jumps
        sys_jumps = ModelingToolkitBase.jumps(compiled_sys)
        @test length(sys_jumps) == 2

        # Both should be VariableRateJumps (state-dependent rates)
        @test all(j -> j isa VariableRateJump, sys_jumps)

        # Poissonians should be cleared after compilation
        @test isempty(ModelingToolkitBase.poissonians(compiled_sys))
    end

    @testset "Constant rate jump" begin
        @parameters λ
        @variables X(t)
        @poissonians dN(λ)

        eqs = [D(X) ~ dN]
        @named sys = System(eqs, t)
        compiled_sys = mtkcompile(sys)

        # Should have 1 jump
        sys_jumps = ModelingToolkitBase.jumps(compiled_sys)
        @test length(sys_jumps) == 1

        # Should be ConstantRateJump (rate is parameter only)
        @test sys_jumps[1] isa ConstantRateJump

        # Equation should be removed (pure-jump)
        @test isempty(equations(compiled_sys))
    end

    @testset "Mixed continuous and jump dynamics" begin
        @parameters a λ
        @variables X(t)
        @poissonians dN(λ)

        # X has continuous decay AND jumps
        eqs = [D(X) ~ -a * X + dN]
        @named sys = System(eqs, t)
        compiled_sys = mtkcompile(sys)

        # Should have 1 jump
        sys_jumps = ModelingToolkitBase.jumps(compiled_sys)
        @test length(sys_jumps) == 1
        @test sys_jumps[1] isa ConstantRateJump

        # Continuous part should be retained
        compiled_eqs = equations(compiled_sys)
        @test length(compiled_eqs) == 1
        # The equation should be D(X) ~ -a*X (poissonian term removed)
    end

    @testset "Time-dependent rate becomes VariableRateJump" begin
        @parameters λ
        @variables X(t)
        @poissonians dN(λ * t)  # Rate depends on t

        eqs = [D(X) ~ dN]
        @named sys = System(eqs, t)
        compiled_sys = mtkcompile(sys)

        sys_jumps = ModelingToolkitBase.jumps(compiled_sys)
        @test length(sys_jumps) == 1
        @test sys_jumps[1] isa VariableRateJump
    end

    @testset "State-dependent rate becomes VariableRateJump" begin
        @parameters k
        @variables X(t)
        @poissonians dN(k * X)  # Rate depends on unknown X

        eqs = [D(X) ~ -dN]
        @named sys = System(eqs, t)
        compiled_sys = mtkcompile(sys)

        sys_jumps = ModelingToolkitBase.jumps(compiled_sys)
        @test length(sys_jumps) == 1
        @test sys_jumps[1] isa VariableRateJump
    end

    @testset "Non-unit jump size (coefficient)" begin
        @parameters λ δ
        @variables X(t)
        @poissonians dN(λ)

        # X jumps by δ, not by 1
        eqs = [D(X) ~ δ * dN]
        @named sys = System(eqs, t)
        compiled_sys = mtkcompile(sys)

        sys_jumps = ModelingToolkitBase.jumps(compiled_sys)
        @test length(sys_jumps) == 1

        # Check that the affect includes the coefficient
        jump = sys_jumps[1]
        @test length(jump.affect!) == 1
    end

    @testset "Same poissonian in multiple equations creates single jump" begin
        @parameters λ
        @variables X(t) Y(t)
        @poissonians dN(λ)

        eqs = [
            D(X) ~ dN,
            D(Y) ~ -dN
        ]
        @named sys = System(eqs, t)
        compiled_sys = mtkcompile(sys)

        # Should create ONE jump with TWO affects
        sys_jumps = ModelingToolkitBase.jumps(compiled_sys)
        @test length(sys_jumps) == 1
        @test length(sys_jumps[1].affect!) == 2
    end

    @testset "Error on nonlinear poissonian usage" begin
        @parameters λ
        @variables X(t)
        @poissonians dN(λ)

        # dN^2 is nonlinear - should error
        eqs = [D(X) ~ dN * dN]
        @named sys = System(eqs, t)

        @test_throws ArgumentError mtkcompile(sys)
    end

    @testset "Merging with explicit jumps" begin
        @parameters λ₁ λ₂
        @variables X(t) Y(t)
        @poissonians dN(λ₁)

        # Create an explicit jump
        explicit_jump = ConstantRateJump(λ₂, [Y ~ Pre(Y) + 1])

        eqs = [D(X) ~ dN]
        @named sys = System(eqs, t; jumps = [explicit_jump])
        compiled_sys = mtkcompile(sys)

        # Should have both jumps
        sys_jumps = ModelingToolkitBase.jumps(compiled_sys)
        @test length(sys_jumps) == 2
    end

    @testset "Various coefficient forms" begin
        @parameters λ a b
        @variables X(t) Y(t) Z(t)
        @poissonians dN(λ)

        # Test dN*a form (poissonian first)
        eqs = [D(X) ~ dN * a]
        @named sys1 = System(eqs, t)
        compiled1 = mtkcompile(sys1)
        @test length(ModelingToolkitBase.jumps(compiled1)) == 1

        # Test (a+b)*dN form (compound coefficient)
        eqs2 = [D(Y) ~ (a + b) * dN]
        @named sys2 = System(eqs2, t)
        compiled2 = mtkcompile(sys2)
        @test length(ModelingToolkitBase.jumps(compiled2)) == 1

        # Test -dN form (negative coefficient)
        eqs3 = [D(Z) ~ -dN]
        @named sys3 = System(eqs3, t)
        compiled3 = mtkcompile(sys3)
        @test length(ModelingToolkitBase.jumps(compiled3)) == 1
    end

    @testset "Jump-diffusion with brownian and poissonian" begin
        @parameters μ σ λ γ
        @variables X(t)
        @brownians dW
        @poissonians dN(λ)

        # dX = μX dt + σX dW + γ dN
        eqs = [D(X) ~ μ * X + σ * X * dW + γ * dN]
        @named sys = System(eqs, t)
        compiled_sys = mtkcompile(sys)

        # Should have 1 jump from poissonian
        sys_jumps = ModelingToolkitBase.jumps(compiled_sys)
        @test length(sys_jumps) == 1
        @test sys_jumps[1] isa ConstantRateJump

        # Should have 1 equation (SDE with drift and diffusion)
        @test length(equations(compiled_sys)) == 1

        # Brownian should be converted to noise equations
        @test !isempty(ModelingToolkitBase.get_noise_eqs(compiled_sys))
    end

    @testset "Hierarchical systems with poissonians" begin
        @parameters λ₁ λ₂
        @variables X(t) Y(t)
        @poissonians dN₁(λ₁) dN₂(λ₂)

        # Create subsystem
        eqs_sub = [D(X) ~ dN₁]
        @named subsys = System(eqs_sub, t)

        # Create parent system
        eqs_parent = [D(Y) ~ dN₂]
        @named parent = System(eqs_parent, t; systems = [subsys])

        # Check poissonians are collected from hierarchy
        all_poissonians = ModelingToolkitBase.poissonians(parent)
        @test length(all_poissonians) == 2

        # Compile and verify both jumps are generated
        compiled = mtkcompile(parent)
        sys_jumps = ModelingToolkitBase.jumps(compiled)
        @test length(sys_jumps) == 2
    end

    @testset "Zero coefficient poissonian is dropped" begin
        @parameters λ₁ λ₂
        @variables X(t)
        @poissonians dN₁(λ₁) dN₂(λ₂)

        # dN₂ has zero coefficient (it appears but with 0 multiplier)
        # Only dN₁ should generate a jump
        eqs = [D(X) ~ dN₁ + 0 * dN₂]
        @named sys = System(eqs, t)
        compiled_sys = mtkcompile(sys)

        # Should only have 1 jump (from dN₁, not dN₂)
        sys_jumps = ModelingToolkitBase.jumps(compiled_sys)
        @test length(sys_jumps) == 1
    end
end
