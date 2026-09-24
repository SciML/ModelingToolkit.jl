using ModelingToolkitBase, Test, JumpProcesses
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using ModelingToolkitBase: ispoissonian, getpoissonianrate, POISSONIAN, getvariabletype
using JumpProcesses: ConstantRateJump, VariableRateJump
using Symbolics: unwrap
using Random, StableRNGs, Statistics
using OrdinaryDiffEq, StochasticDiffEq

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

    @testset "Rate from local variable expression" begin
        # Test that a locally computed rate expression works
        @parameters a b
        rate_expr = a * b^2
        @poissonians dN(rate_expr)

        # Rate should be the evaluated symbolic expression, not the variable name
        @test isequal(getpoissonianrate(dN), a * b^2)
    end

    @testset "Non-poissonian returns nothing for rate" begin
        @variables x(t)

        # Regular variables should return nothing for rate
        @test getpoissonianrate(x) === nothing
        @test !ispoissonian(x)
    end

    @testset "Return type consistency with @variables" begin
        @parameters λ₁ λ₂ λ₃

        # Single inline: should return a 1-element Vector
        result_single = @poissonians dN_a(λ₁)
        @test result_single isa Vector
        @test length(result_single) == 1
        @test ispoissonian(result_single[1])
        @test isequal(result_single[1], dN_a)

        # Multiple inline: should return a Vector
        result_multi = @poissonians dN_b(λ₁) dN_c(λ₂)
        @test result_multi isa Vector
        @test length(result_multi) == 2
        @test all(ispoissonian, result_multi)
        @test isequal(result_multi[1], dN_b)
        @test isequal(result_multi[2], dN_c)

        # Single in block: should return a 1-element Vector
        result_block_single = @poissonians begin
            dN_d(λ₁)
        end
        @test result_block_single isa Vector
        @test length(result_block_single) == 1
        @test ispoissonian(result_block_single[1])
        @test isequal(result_block_single[1], dN_d)

        # Multiple in block: should return a Vector
        result_block_multi = @poissonians begin
            dN_e(λ₂)
            dN_f(λ₃)
        end
        @test result_block_multi isa Vector
        @test length(result_block_multi) == 2
        @test all(ispoissonian, result_block_multi)
        @test isequal(result_block_multi[1], dN_e)
        @test isequal(result_block_multi[2], dN_f)

        # Three inline: should return a Vector of length 3
        result_three = @poissonians dN_g(λ₁) dN_h(λ₂) dN_i(λ₃)
        @test result_three isa Vector
        @test length(result_three) == 3
        @test all(ispoissonian, result_three)
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
            D(R) ~ dN_rec,
        ]

        # Using explicit 5-argument constructor with poissonians kwarg
        @named sys = System(eqs, t, [S, I, R], [β, γ]; poissonians = [dN_inf, dN_rec])

        @test length(ModelingToolkitBase.poissonians(sys)) == 2
        @test Set(ModelingToolkitBase.poissonians(sys)) == Set((dN_inf, dN_rec))
    end

    @testset "Two-argument constructor auto-detects poissonians" begin
        @parameters β γ
        @variables S(t) I(t) R(t)
        @poissonians dN_inf(β * S * I) dN_rec(γ * I)

        eqs = [
            D(S) ~ -dN_inf,
            D(I) ~ dN_inf - dN_rec,
            D(R) ~ dN_rec,
        ]

        # Two-argument constructor should auto-detect poissonians
        @named sys = System(eqs, t)

        @test length(ModelingToolkitBase.poissonians(sys)) == 2
        @test Set(ModelingToolkitBase.poissonians(sys)) == Set((dN_inf, dN_rec))
    end

    @testset "Variables from rate expressions are extracted" begin
        @parameters β γ
        @variables S(t) I(t) R(t)
        @poissonians dN_inf(β * S * I) dN_rec(γ * I)

        eqs = [
            D(S) ~ -dN_inf,
            D(I) ~ dN_inf - dN_rec,
            D(R) ~ dN_rec,
        ]

        @named sys = System(eqs, t)

        # Exactly S, I, R should be in unknowns (no poissonians or extras)
        @test Set(ModelingToolkitBase.unknowns(sys)) == Set((S, I, R))

        # Exactly β, γ should be in parameters (from rate expressions)
        @test Set(ModelingToolkitBase.parameters(sys)) == Set((β, γ))
    end

    @testset "Unknown in coefficient is auto-discovered" begin
        @parameters λ
        @variables X(t)
        @poissonians dN(λ)

        # X appears in coefficient, not in rate expression
        eqs = [D(X) ~ X * dN]

        @named sys = System(eqs, t)

        # Exactly X should be in unknowns (no poissonian)
        @test Set(ModelingToolkitBase.unknowns(sys)) == Set((X,))

        # Exactly λ should be in parameters (from rate expression)
        @test Set(ModelingToolkitBase.parameters(sys)) == Set((λ,))
    end

    @testset "Parameter in coefficient is auto-discovered" begin
        @parameters λ α
        @variables X(t)
        @poissonians dN(λ)

        # α appears in coefficient, λ in rate
        eqs = [D(X) ~ α * dN]

        @named sys = System(eqs, t)

        # Exactly λ and α should be in parameters
        @test Set(ModelingToolkitBase.parameters(sys)) == Set((λ, α))

        # Exactly X should be in unknowns (from LHS, no poissonian)
        @test Set(ModelingToolkitBase.unknowns(sys)) == Set((X,))
    end

    @testset "Mixed parameter and unknown in coefficient auto-discovered" begin
        @parameters λ α β
        @variables X(t) Y(t)
        @poissonians dN(λ)

        # α*X + β*Y in coefficient
        eqs = [D(X) ~ (α * X + β * Y) * dN]

        @named sys = System(eqs, t)

        # Exactly λ, α, β should be in parameters
        @test Set(ModelingToolkitBase.parameters(sys)) == Set((λ, α, β))

        # Exactly X and Y should be in unknowns (no poissonian)
        @test Set(ModelingToolkitBase.unknowns(sys)) == Set((X, Y))
    end

    @testset "Poissonians are not included in unknowns" begin
        @parameters λ
        @variables X(t)
        @poissonians dN(λ)

        eqs = [D(X) ~ dN]

        @named sys = System(eqs, t)

        # Exactly X should be in unknowns (dN should not be)
        @test Set(ModelingToolkitBase.unknowns(sys)) == Set((X,))
        # Exactly dN should be in poissonians
        @test Set(ModelingToolkitBase.poissonians(sys)) == Set((dN,))
        # Exactly λ should be in parameters
        @test Set(ModelingToolkitBase.parameters(sys)) == Set((λ,))
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
            D(R) ~ dN_rec,
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
            D(Y) ~ -dN,
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

    @testset "Error on higher-order derivatives with poissonians" begin
        @parameters λ
        @variables X(t) Y(t)
        @poissonians dN(λ)

        # Any higher-order derivative in a system with poissonians should error
        # (model is not well-defined with jumps and higher-order derivatives)
        eqs = [D(D(X)) ~ dN]
        @named sys = System(eqs, t)
        @test_throws ArgumentError mtkcompile(sys)

        # Even if poissonian is in a different equation, higher-order derivative errors
        eqs2 = [D(D(X)) ~ 0, D(Y) ~ dN]
        @named sys2 = System(eqs2, t)
        @test_throws ArgumentError mtkcompile(sys2)
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

    @testset "State-dependent coefficient (unknown in coefficient)" begin
        @parameters λ
        @variables X(t)
        @poissonians dN(λ)

        # Jump size is proportional to current state: X jumps by X (i.e., doubles)
        eqs = [D(X) ~ X * dN]
        @named sys = System(eqs, t)
        compiled_sys = mtkcompile(sys)

        sys_jumps = ModelingToolkitBase.jumps(compiled_sys)
        @test length(sys_jumps) == 1
        @test sys_jumps[1] isa ConstantRateJump  # Rate is λ (parameter only)

        # Verify the affect exists and targets X
        @test length(sys_jumps[1].affect!) == 1
        affect_eq = sys_jumps[1].affect![1]
        @test isequal(affect_eq.lhs, X)
    end

    @testset "Mixed parameter and unknown coefficient" begin
        @parameters λ α
        @variables X(t)
        @poissonians dN(λ)

        # Jump size is α*X (parameter times state)
        eqs = [D(X) ~ α * X * dN]
        @named sys = System(eqs, t)
        compiled_sys = mtkcompile(sys)

        sys_jumps = ModelingToolkitBase.jumps(compiled_sys)
        @test length(sys_jumps) == 1

        # Verify the affect exists and targets X
        @test length(sys_jumps[1].affect!) == 1
        affect_eq = sys_jumps[1].affect![1]
        @test isequal(affect_eq.lhs, X)
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

    @testset "save_positions forwarded to VRJs from poissonians" begin
        @parameters k
        @variables X(t)
        @poissonians dN(k * X)  # State-dependent → VRJ

        eqs = [D(X) ~ -dN]
        @named sys = System(eqs, t)

        # Default save_positions = (false, true)
        compiled_default = mtkcompile(sys)
        vrj_default = ModelingToolkitBase.jumps(compiled_default)[1]
        @test vrj_default isa VariableRateJump
        @test vrj_default.save_positions == (false, true)

        # Custom save_positions = (true, false)
        compiled_custom = mtkcompile(sys; save_positions = (true, false))
        vrj_custom = ModelingToolkitBase.jumps(compiled_custom)[1]
        @test vrj_custom isa VariableRateJump
        @test vrj_custom.save_positions == (true, false)

        # Custom save_positions = (false, false)
        compiled_none = mtkcompile(sys; save_positions = (false, false))
        vrj_none = ModelingToolkitBase.jumps(compiled_none)[1]
        @test vrj_none.save_positions == (false, false)
    end

    @testset "save_positions does not affect user-provided VRJs" begin
        @parameters λ k
        @variables X(t) Y(t)
        @poissonians dN(λ)  # Constant rate → CRJ (no save_positions field)

        # User-provided VRJ with explicit save_positions
        user_vrj = VariableRateJump(k * Y, [Y ~ Pre(Y) - 1]; save_positions = (true, true))

        eqs = [D(X) ~ dN]
        @named sys = System(eqs, t; jumps = [user_vrj])

        # Compile with different save_positions - should NOT affect user VRJ
        compiled = mtkcompile(sys; save_positions = (false, false))
        sys_jumps = ModelingToolkitBase.jumps(compiled)

        # Find the user VRJ (it's the one with rate k*Y)
        user_vrj_compiled = filter(j -> j isa VariableRateJump, sys_jumps)
        @test length(user_vrj_compiled) == 1
        @test user_vrj_compiled[1].save_positions == (true, true)  # Unchanged!
    end

    @testset "Error when JumpProblem receives uncompiled system with poissonians" begin
        @parameters λ
        @variables X(t)
        @poissonians dN(λ)

        eqs = [D(X) ~ dN]
        @named sys = System(eqs, t)

        # complete() but NOT mtkcompile() - poissonians are still unprocessed
        sys_complete = complete(sys)

        # Should error because poissonians haven't been converted to jumps
        @test_throws ModelingToolkitBase.SystemCompatibilityError JumpProblem(
            sys_complete, [X => 0.0, λ => 1.0], (0.0, 1.0)
        )
    end
end

# ============================================================================
# Phase 4: Mathematical Correctness / Integration Tests
# ============================================================================

@testset "Poissonian Mathematical Correctness" begin
    rng = StableRNG(12345)

    @testset "Pure Poisson counter: E[N(t)] = λt, Var[N(t)] = λt" begin
        # Simple counting process: dN = dN(λ)
        # N(t) ~ Poisson(λt), so E[N(t)] = λt, Var[N(t)] = λt
        # Rate is parameter-only → ConstantRateJump → can use SSAStepper
        @parameters λ
        @variables N(t)
        @poissonians dN(λ)

        eqs = [D(N) ~ dN]
        @named sys = System(eqs, t)
        compiled_sys = mtkcompile(sys)

        # Verify we got a ConstantRateJump (rate depends only on parameters)
        sys_jumps = ModelingToolkitBase.jumps(compiled_sys)
        @test length(sys_jumps) == 1
        @test sys_jumps[1] isa ConstantRateJump

        λ_val = 5.0
        T = 2.0
        Nsims = 2000

        # Pure jump system with ConstantRateJump → use SSAStepper
        # SSAStepper accepts seed kwarg for reproducibility (not save_everystep!)
        jprob = JumpProblem(
            compiled_sys, [N => 0.0, λ => λ_val], (0.0, T);
            aggregator = Direct(), save_positions = (false, false), rng
        )

        seed = 1111
        Nfinal = zeros(Nsims)
        for i in 1:Nsims
            sol = solve(jprob, SSAStepper(); seed)
            Nfinal[i] = sol[N, end]
            seed += 1
        end

        E_N = λ_val * T  # = 10
        Var_N = λ_val * T  # = 10

        sample_mean = mean(Nfinal)
        sample_var = var(Nfinal)

        @test abs(sample_mean - E_N) < 0.05 * E_N  # 5% relative error
        @test abs(sample_var - Var_N) < 0.1 * Var_N  # 10% for variance
    end

    @testset "Birth-death process: compare to analytical steady state" begin
        # Birth at constant rate b, death at state-dependent rate d*X
        # State-dependent rate → VariableRateJump → needs ODE solver (Tsit5)
        # E[X(∞)] = b/d (steady state)
        @parameters b d
        @variables X(t)
        @poissonians dN_birth(b) dN_death(d * X)

        eqs = [D(X) ~ dN_birth - dN_death]
        @named sys = System(eqs, t)
        compiled_sys = mtkcompile(sys)

        # Verify we have both CRJ and VRJ (birth is constant, death is variable)
        sys_jumps = ModelingToolkitBase.jumps(compiled_sys)
        @test length(sys_jumps) == 2
        @test count(j -> j isa ConstantRateJump, sys_jumps) == 1
        @test count(j -> j isa VariableRateJump, sys_jumps) == 1

        b_val = 5.0
        d_val = 0.1
        X0 = 10.0
        T = 100.0  # Long time to reach steady state
        Nsims = 1000

        # VariableRateJumps require ODE solver
        jprob_sym = JumpProblem(
            compiled_sys, [X => X0, b => b_val, d => d_val], (0.0, T);
            save_positions = (false, false), rng
        )

        seed = 2222
        Xfinal_sym = zeros(Nsims)
        for i in 1:Nsims
            sol = solve(jprob_sym, Tsit5(); save_everystep = false, seed)
            Xfinal_sym[i] = sol[X, end]
            seed += 1
        end

        E_X_ss = b_val / d_val  # = 50

        mean_sym = mean(Xfinal_sym)

        # Should match analytical steady state
        @test abs(mean_sym - E_X_ss) < 0.1 * E_X_ss
    end

    @testset "SIR model: compare symbolic @poissonians to direct JumpProcesses" begin
        # Classic SIR: S + I --(β*S*I)--> 2I, I --(γ*I)--> R
        # State-dependent rates → VariableRateJumps → needs ODE solver
        @parameters β γ
        @variables S(t) I(t) R(t)
        @poissonians dN_inf(β * S * I) dN_rec(γ * I)

        eqs = [
            D(S) ~ -dN_inf,
            D(I) ~ dN_inf - dN_rec,
            D(R) ~ dN_rec,
        ]
        @named sys = System(eqs, t)
        compiled_sys = mtkcompile(sys)

        # Both rates depend on unknowns → both are VariableRateJumps
        sys_jumps = ModelingToolkitBase.jumps(compiled_sys)
        @test length(sys_jumps) == 2
        @test all(j -> j isa VariableRateJump, sys_jumps)

        β_val = 0.1 / 1000
        γ_val = 0.01
        S0, I0, R0 = 999.0, 1.0, 0.0
        T = 250.0
        Nsims = 2000
        save_times = 1.0:1.0:T
        Ntimes = length(save_times)

        jprob_sym = JumpProblem(
            compiled_sys,
            [S => S0, I => I0, R => R0, β => β_val, γ => γ_val], (0.0, T);
            save_positions = (false, false), rng
        )

        seed = 3333
        R_sym = zeros(Nsims, Ntimes)
        for i in 1:Nsims
            sol = solve(jprob_sym, Tsit5(); seed)
            for (k, t_) in enumerate(save_times)
                R_sym[i, k] = sol(t_; idxs = R)
            end
            seed += 1
        end

        # Build directly with JumpProcesses using VariableRateJumps for fair comparison
        f_direct(du, u, p, t) = (du .= 0)  # No continuous dynamics
        oprob_direct = ODEProblem(f_direct, [S0, I0, R0], (0.0, T), (β_val, γ_val))
        r1(u, p, t) = p[1] * u[1] * u[2]
        function a1!(integ)
            integ.u[1] -= 1
            integ.u[2] += 1
        end
        j1 = VariableRateJump(r1, a1!)
        r2(u, p, t) = p[2] * u[2]
        function a2!(integ)
            integ.u[2] -= 1
            integ.u[3] += 1
        end
        j2 = VariableRateJump(r2, a2!)
        jprob_direct = JumpProblem(
            oprob_direct, Direct(), j1, j2;
            rng, save_positions = (false, false)
        )

        seed = 3333
        R_direct = zeros(Nsims, Ntimes)
        for i in 1:Nsims
            sol = solve(jprob_direct, Tsit5(); seed)
            for (k, t_) in enumerate(save_times)
                R_direct[i, k] = sol(t_)[3]
            end
            seed += 1
        end

        # Compare means at every saved time point
        mean_sym = vec(mean(R_sym; dims = 1))
        mean_direct = vec(mean(R_direct; dims = 1))

        # All time points should match within tolerance
        @test all(
            abs(mean_sym[k] - mean_direct[k]) < 0.1 * max(mean_direct[k], 1.0)
                for k in 1:Ntimes
        )
    end

    @testset "Jump-diffusion: @brownians + @poissonians vs analytical" begin
        # dX = σ*dW + δ*dN(λ)
        # E[X(T)] = δ*λ*T (diffusion has zero mean)
        # Var[X(T)] = σ²*T + δ²*λ*T
        @parameters σ λ δ
        @variables X(t)
        @brownians dW
        @poissonians dN(λ)

        eqs = [D(X) ~ σ * dW + δ * dN]
        @named sys = System(eqs, t)
        compiled_sys = mtkcompile(sys)

        σ_val = 0.3
        λ_val = 2.0
        δ_val = 1.0
        T = 2.0
        Nsims = 2000

        jprob_sym = JumpProblem(
            compiled_sys,
            [X => 0.0, σ => σ_val, λ => λ_val, δ => δ_val], (0.0, T);
            save_positions = (false, false), rng
        )

        seed = 4444
        Xfinal_sym = zeros(Nsims)
        for i in 1:Nsims
            sol = solve(jprob_sym, SOSRI(); save_everystep = false, seed)
            Xfinal_sym[i] = sol[X, end]
            seed += 1
        end

        E_X = δ_val * λ_val * T  # = 4.0
        Var_X = σ_val^2 * T + δ_val^2 * λ_val * T  # = 0.18 + 4 = 4.18

        sample_mean = mean(Xfinal_sym)
        sample_var = var(Xfinal_sym)

        @test abs(sample_mean - E_X) < 0.1 * E_X
        @test abs(sample_var - Var_X) < 0.15 * Var_X
    end

    @testset "Mixed continuous + jump: decay with immigration" begin
        # dX = -a*X*dt + dN(λ)
        # At steady state: E[X(∞)] = λ/a
        # λ is parameter-only → ConstantRateJump, but has ODE part
        @parameters a λ
        @variables X(t)
        @poissonians dN(λ)

        eqs = [D(X) ~ -a * X + dN]
        @named sys = System(eqs, t)
        compiled_sys = mtkcompile(sys)

        # Jump should be ConstantRateJump (rate is just λ)
        sys_jumps = ModelingToolkitBase.jumps(compiled_sys)
        @test length(sys_jumps) == 1
        @test sys_jumps[1] isa ConstantRateJump

        a_val = 0.5
        λ_val = 5.0
        X0 = 0.0
        T = 50.0  # Long time to reach steady state
        Nsims = 1000

        # Has ODE part, so needs ODE solver
        jprob = JumpProblem(
            compiled_sys, [X => X0, a => a_val, λ => λ_val], (0.0, T);
            save_positions = (false, false), rng
        )

        seed = 5555
        Xfinal = zeros(Nsims)
        for i in 1:Nsims
            sol = solve(jprob, Tsit5(); save_everystep = false, seed)
            Xfinal[i] = sol[X, end]
            seed += 1
        end

        E_X_ss = λ_val / a_val  # = 10

        sample_mean = mean(Xfinal)
        @test abs(sample_mean - E_X_ss) < 0.1 * E_X_ss
    end

    @testset "State-dependent coefficient: geometric decay X ~ -δ*X*dN(λ)" begin
        # Jump size is -δ*X (proportional decay)
        # After n jumps: X(t) = X₀ * (1-δ)^n
        # E[X(T)] = X₀ * E[(1-δ)^N] where N ~ Poisson(λT)
        # Using MGF of Poisson: E[e^{θN}] = exp(λT(e^θ - 1))
        # Set θ = ln(1-δ): E[(1-δ)^N] = exp(λT((1-δ) - 1)) = exp(-λTδ)
        # So: E[X(T)] = X₀ * exp(-λ*δ*T)
        @parameters λ δ
        @variables X(t)
        @poissonians dN(λ)

        # Jump decreases X by δ*X (fractional decay)
        eqs = [D(X) ~ -δ * X * dN]
        @named sys = System(eqs, t)
        compiled_sys = mtkcompile(sys)

        # Jump should be ConstantRateJump (rate is λ, parameter-only)
        # but the affect has state-dependent coefficient
        sys_jumps = ModelingToolkitBase.jumps(compiled_sys)
        @test length(sys_jumps) == 1
        @test sys_jumps[1] isa ConstantRateJump

        λ_val = 2.0
        δ_val = 0.2  # 20% decay per jump
        X0 = 100.0
        T = 3.0
        Nsims = 2000

        # Pure jump system with CRJ → can use SSAStepper
        jprob = JumpProblem(
            compiled_sys, [X => X0, λ => λ_val, δ => δ_val], (0.0, T);
            aggregator = Direct(), save_positions = (false, false), rng
        )

        seed = 6666
        Xfinal = zeros(Nsims)
        for i in 1:Nsims
            sol = solve(jprob, SSAStepper(); seed)
            Xfinal[i] = sol[X, end]
            seed += 1
        end

        # Analytical: E[X(T)] = X₀ * exp(-λ*δ*T)
        E_X_T = X0 * exp(-λ_val * δ_val * T)  # = 100 * exp(-1.2) ≈ 30.12

        sample_mean = mean(Xfinal)

        @test abs(sample_mean - E_X_T) < 0.1 * E_X_T
    end

    @testset "State-dependent coefficient: compare to direct JumpProcesses" begin
        # Same geometric decay model, compare symbolic to direct JumpProcesses
        @parameters λ δ
        @variables X(t)
        @poissonians dN(λ)

        eqs = [D(X) ~ -δ * X * dN]
        @named sys = System(eqs, t)
        compiled_sys = mtkcompile(sys)

        λ_val = 3.0
        δ_val = 0.15
        X0 = 50.0
        T = 2.0
        Nsims = 1000

        # Symbolic version
        jprob_sym = JumpProblem(
            compiled_sys, [X => X0, λ => λ_val, δ => δ_val], (0.0, T);
            aggregator = Direct(), save_positions = (false, false), rng
        )

        seed = 7777
        Xfinal_sym = zeros(Nsims)
        for i in 1:Nsims
            sol = solve(jprob_sym, SSAStepper(); seed)
            Xfinal_sym[i] = sol[X, end]
            seed += 1
        end

        # Direct JumpProcesses version
        u0 = [X0]
        dprob = DiscreteProblem(u0, (0.0, T), (λ_val, δ_val))
        rate_direct(u, p, t) = p[1]  # λ
        function affect_direct!(integ)
            integ.u[1] -= integ.p[2] * integ.u[1]  # X -= δ*X
        end
        crj_direct = ConstantRateJump(rate_direct, affect_direct!)
        jprob_direct = JumpProblem(dprob, Direct(), crj_direct; rng)

        seed = 7777
        Xfinal_direct = zeros(Nsims)
        for i in 1:Nsims
            sol = solve(jprob_direct, SSAStepper(); seed)
            Xfinal_direct[i] = sol[end][1]
            seed += 1
        end

        mean_sym = mean(Xfinal_sym)
        mean_direct = mean(Xfinal_direct)

        # Both should give same results (same seeds)
        @test abs(mean_sym - mean_direct) < 0.05 * mean_direct
    end
end
