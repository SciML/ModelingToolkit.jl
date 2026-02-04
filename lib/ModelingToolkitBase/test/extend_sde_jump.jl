using ModelingToolkitBase, Test
using StochasticDiffEq, JumpProcesses, OrdinaryDiffEq
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using Statistics
using StableRNGs
using Symbolics: value

MT = ModelingToolkitBase
rng = StableRNG(12345)

@testset "Extend with brownians" begin
    @testset "Extend two systems with disjoint brownians" begin
        @variables x(t) y(t)
        @parameters a b
        @brownians w1 w2

        @named sys1 = System([D(x) ~ a * x + w1 * x], t, [x], [a], [w1])
        @named sys2 = System([D(y) ~ b * y + w2 * y], t, [y], [b], [w2])

        extended = extend(sys2, sys1; name = :extended)

        # Verify brownians are merged
        @test length(MT.get_brownians(extended)) == 2
        @test Set(MT.get_brownians(extended)) == Set([w1, w2])

        # Verify can compile and solve
        compiled_sys = mtkcompile(extended)
        prob = SDEProblem(compiled_sys, [x => 1.0, y => 1.0, a => 0.1, b => 0.2], (0.0, 1.0))
        sol = solve(prob, EM(); dt = 0.01)
        @test sol.retcode == ReturnCode.Success
    end

    @testset "Extend with shared brownian" begin
        @variables x(t) y(t)
        @parameters a b
        @brownians w

        @named sys1 = System([D(x) ~ a * x + w * x], t, [x], [a], [w])
        @named sys2 = System([D(y) ~ b * y + w * y], t, [y], [b], [w])

        extended = extend(sys2, sys1; name = :extended)

        # Verify brownians are unioned (no duplicates)
        @test length(MT.get_brownians(extended)) == 1
        @test isequal(MT.get_brownians(extended)[1], w)

        # Verify can compile and solve - both variables should receive same kicks
        compiled_sys = mtkcompile(extended)
        prob = SDEProblem(compiled_sys, [x => 1.0, y => 1.0, a => 0.0, b => 0.0], (0.0, 1.0))
        sol = solve(prob, EM(); dt = 0.01)
        @test sol.retcode == ReturnCode.Success
        # With zero drift and same brownian, final values should be equal
        @test sol[x][end] == sol[y][end]
    end

    @testset "Extend base ODE with SDE system" begin
        @variables x(t) y(t)
        @parameters a b
        @brownians w

        @named base_ode = System([D(x) ~ a * x], t, [x], [a])
        @named sde_sys = System([D(y) ~ b * y + w * y], t, [y], [b], [w])

        extended = extend(sde_sys, base_ode; name = :extended)

        # Verify brownians from SDE system are present
        @test length(MT.get_brownians(extended)) == 1

        # Verify can compile as SDE
        compiled_sys = mtkcompile(extended)
        prob = SDEProblem(compiled_sys, [x => 1.0, y => 1.0, a => 0.1, b => 0.2], (0.0, 1.0))
        sol = solve(prob, EM(); dt = 0.01)
        @test sol.retcode == ReturnCode.Success
    end
end

@testset "Extend with poissonians" begin
    @testset "Extend two systems with disjoint poissonians (constant rate)" begin
        @variables X(t) Y(t)
        @parameters λ1 λ2
        @poissonians dN1(λ1) dN2(λ2)

        @named sys1 = System([D(X) ~ dN1], t, [X], [λ1]; poissonians = [dN1])
        @named sys2 = System([D(Y) ~ dN2], t, [Y], [λ2]; poissonians = [dN2])

        extended = extend(sys2, sys1; name = :extended)

        # Verify poissonians are merged
        @test length(MT.get_poissonians(extended)) == 2

        # Verify can compile - poissonians become jumps
        compiled_sys = mtkcompile(extended)
        @test length(MT.jumps(compiled_sys)) == 2

        # Verify can solve
        jprob = JumpProblem(compiled_sys, [X => 0, Y => 0, λ1 => 1.0, λ2 => 2.0], (0.0, 10.0); rng)
        sol = solve(jprob, SSAStepper())
        @test sol.retcode == ReturnCode.Success
    end

    @testset "Extend with shared poissonian" begin
        @variables X(t) Y(t)
        @parameters λ
        @poissonians dN(λ)

        @named sys1 = System([D(X) ~ dN], t, [X], [λ]; poissonians = [dN])
        @named sys2 = System([D(Y) ~ 2 * dN], t, [Y], []; poissonians = [dN])

        extended = extend(sys2, sys1; name = :extended)

        # Verify poissonians are unioned (no duplicates)
        @test length(MT.get_poissonians(extended)) == 1

        # Verify can compile - should produce single jump with both affects
        compiled_sys = mtkcompile(extended)
        @test length(MT.jumps(compiled_sys)) == 1

        # Verify can solve
        jprob = JumpProblem(compiled_sys, [X => 0, Y => 0, λ => 2.0], (0.0, 10.0); rng)
        sol = solve(jprob, SSAStepper())
        @test sol.retcode == ReturnCode.Success
        # Y should increase twice as fast as X
        @test sol[Y][end] ≈ 2 * sol[X][end]
    end

    @testset "Extend ODE with poissonian system" begin
        @variables X(t) Y(t)
        @parameters a λ
        @poissonians dN(λ)

        @named base_ode = System([D(X) ~ a * X], t, [X], [a])
        @named pois_sys = System([D(Y) ~ dN], t, [Y], [λ]; poissonians = [dN])

        extended = extend(pois_sys, base_ode; name = :extended)

        @test length(MT.get_poissonians(extended)) == 1
        @test length(MT.get_eqs(extended)) == 2

        # Verify can compile - produces hybrid ODE + jump system
        compiled_sys = mtkcompile(extended)
        @test length(MT.jumps(compiled_sys)) == 1
    end

    @testset "Extend with variable rate poissonian" begin
        @variables X(t) Y(t)
        @parameters k
        @poissonians dN(k * X)  # Rate depends on unknown -> VariableRateJump

        # ODE part with the unknown X that the poissonian rate depends on
        @named sys1 = System([D(X) ~ -0.1 * X], t, [X], [k])
        # Poissonian part - equation using dN so a jump is generated
        @named sys2 = System([D(Y) ~ dN], t, [Y], []; poissonians = [dN])

        extended = extend(sys2, sys1; name = :extended)

        @test length(MT.get_poissonians(extended)) == 1

        # Verify can compile - should produce VariableRateJump
        compiled_sys = mtkcompile(extended)
        @test length(MT.jumps(compiled_sys)) == 1
        @test MT.jumps(compiled_sys)[1] isa VariableRateJump
    end
end

@testset "Extend with jumps" begin
    @testset "Extend two systems with disjoint ConstantRateJumps" begin
        @variables S(t) I(t) R(t)
        @parameters β γ

        # First jump system: S -> I
        rate1 = β * S * I
        affect1 = [S ~ Pre(S) - 1, I ~ Pre(I) + 1]
        j1 = ConstantRateJump(rate1, affect1)
        @named sys1 = JumpSystem([j1], t, [S, I], [β])

        # Second jump system: I -> R
        rate2 = γ * I
        affect2 = [I ~ Pre(I) - 1, R ~ Pre(R) + 1]
        j2 = ConstantRateJump(rate2, affect2)
        @named sys2 = JumpSystem([j2], t, [I, R], [γ])

        extended = extend(sys2, sys1; name = :sir)

        # Verify jumps are merged
        @test length(MT.get_jumps(extended)) == 2

        # Verify can solve
        completed = complete(extended)
        u0 = [S => 99, I => 1, R => 0]
        ps = [β => 0.1 / 100, γ => 0.05]
        jprob = JumpProblem(completed, [u0; ps], (0.0, 100.0); rng)
        sol = solve(jprob, SSAStepper())
        @test sol.retcode == ReturnCode.Success
        # Conservation check: S + I + R should equal 100
        @test sol[S][end] + sol[I][end] + sol[R][end] == 100
    end

    @testset "Extend with MassActionJumps" begin
        @variables A(t) B(t)
        @parameters k1 k2

        maj1 = MassActionJump(k1, [A => 1], [A => -1, B => 1])
        @named sys1 = JumpSystem([maj1], t, [A, B], [k1])

        maj2 = MassActionJump(k2, [B => 1], [B => -1, A => 1])
        @named sys2 = JumpSystem([maj2], t, [A, B], [k2])

        extended = extend(sys2, sys1; name = :extended)

        @test length(MT.get_jumps(extended)) == 2

        completed = complete(extended)
        u0 = [A => 50, B => 50]
        ps = [k1 => 0.1, k2 => 0.1]
        jprob = JumpProblem(completed, [u0; ps], (0.0, 10.0); rng)
        sol = solve(jprob, SSAStepper())
        @test sol.retcode == ReturnCode.Success
        # Conservation: A + B = 100
        @test sol[A][end] + sol[B][end] == 100
    end

    @testset "Extend ODE base with jump system" begin
        @variables X(t) S(t)
        @parameters a β

        @named base_ode = System([D(X) ~ a * X], t, [X], [a])

        rate = β * S
        affect = [S ~ Pre(S) - 1]
        j = ConstantRateJump(rate, affect)
        @named jump_sys = JumpSystem([j], t, [S], [β])

        extended = extend(jump_sys, base_ode; name = :extended)

        @test length(MT.get_jumps(extended)) == 1
        @test length(MT.get_eqs(extended)) == 1  # ODE equation
    end
end

@testset "Extend with mixed stochastic components" begin
    @testset "Extend brownians + poissonians" begin
        @variables X(t) Y(t)
        @parameters σ λ
        @brownians w
        @poissonians dN(λ)

        @named sde_sys = System([D(X) ~ σ * w], t, [X], [σ], [w])
        @named pois_sys = System([D(Y) ~ dN], t, [Y], [λ]; poissonians = [dN])

        extended = extend(pois_sys, sde_sys; name = :extended)

        @test length(MT.get_brownians(extended)) == 1
        @test length(MT.get_poissonians(extended)) == 1

        # Verify can compile - produces SDE with jumps
        compiled_sys = mtkcompile(extended)
        @test MT.get_noise_eqs(compiled_sys) !== nothing
        @test length(MT.jumps(compiled_sys)) == 1
    end

    @testset "Extend brownians + jumps" begin
        @variables X(t) Y(t)
        @parameters σ k
        @brownians w

        @named sde_sys = System([D(X) ~ σ * w], t, [X], [σ], [w])

        rate = k
        affect = [Y ~ Pre(Y) + 1]
        j = ConstantRateJump(rate, affect)
        @named jump_sys = JumpSystem([j], t, [Y], [k])

        extended = extend(jump_sys, sde_sys; name = :extended)

        @test length(MT.get_brownians(extended)) == 1
        @test length(MT.get_jumps(extended)) == 1
    end

    @testset "Extend poissonians + jumps" begin
        @variables X(t) Y(t)
        @parameters λ k
        @poissonians dN(λ)

        @named pois_sys = System([D(X) ~ dN], t, [X], [λ]; poissonians = [dN])

        rate = k
        affect = [Y ~ Pre(Y) + 1]
        j = ConstantRateJump(rate, affect)
        @named jump_sys = JumpSystem([j], t, [Y], [k])

        extended = extend(jump_sys, pois_sys; name = :extended)

        @test length(MT.get_poissonians(extended)) == 1
        @test length(MT.get_jumps(extended)) == 1

        # After compile, should have 2 jumps (1 from poissonian + 1 explicit)
        compiled_sys = mtkcompile(extended)
        @test length(MT.jumps(compiled_sys)) == 2
    end

    @testset "Extend all three: brownians + poissonians + jumps" begin
        @variables X(t) Y(t) Z(t)
        @parameters σ λ k
        @brownians w
        @poissonians dN(λ)

        @named sde_sys = System([D(X) ~ σ * w], t, [X], [σ], [w])
        @named pois_sys = System([D(Y) ~ dN], t, [Y], [λ]; poissonians = [dN])

        rate = k
        affect = [Z ~ Pre(Z) + 1]
        j = ConstantRateJump(rate, affect)
        @named jump_sys = JumpSystem([j], t, [Z], [k])

        # Chain extend: sde + pois + jump
        extended1 = extend(pois_sys, sde_sys; name = :temp)
        extended = extend(jump_sys, extended1; name = :all_three)

        @test length(MT.get_brownians(extended)) == 1
        @test length(MT.get_poissonians(extended)) == 1
        @test length(MT.get_jumps(extended)) == 1
    end
end

@testset "Priority semantics" begin
    @testset "Initial conditions prefer sys over basesys" begin
        @variables x(t)
        @parameters a
        @brownians w

        # Base system with initial condition
        @named base = System([D(x) ~ a * x + w], t, [x], [a], [w];
            initial_conditions = Dict(x => 1.0))

        # Extending system with different initial condition
        @named ext = System(Equation[], t, [], [], [];
            initial_conditions = Dict(x => 2.0))

        extended = extend(ext, base; name = :result)

        # Initial condition should prefer ext's value (2.0)
        @test value(MT.get_initial_conditions(extended)[x]) == 2.0
    end

    @testset "Brownians from both systems are unioned" begin
        @variables x(t) y(t)
        @brownians w1 w2

        @named sys1 = System([D(x) ~ w1], t, [x], [], [w1])
        @named sys2 = System([D(y) ~ w2], t, [y], [], [w2])

        extended = extend(sys2, sys1; name = :result)

        # Both brownians should be present
        @test length(MT.get_brownians(extended)) == 2
        @test Set(MT.get_brownians(extended)) == Set([w1, w2])
    end
end

@testset "Correctness tests" begin
    @testset "Extended SDE statistics - shared brownian" begin
        # Two variables with same brownian should have correlated paths
        @variables x(t) y(t)
        @brownians w

        @named sys1 = System([D(x) ~ w], t, [x], [], [w])
        @named sys2 = System([D(y) ~ w], t, [y], [], [w])

        extended = extend(sys2, sys1; name = :extended)
        compiled_sys = mtkcompile(extended)

        # Run many simulations
        N = 100
        equal_count = 0
        for i in 1:N
            prob = SDEProblem(compiled_sys, [x => 0.0, y => 0.0], (0.0, 1.0))
            sol = solve(prob, EM(); dt = 0.001, save_everystep = false)
            # Since both use same brownian with coefficient 1, they should be equal
            if sol[x][end] ≈ sol[y][end]
                equal_count += 1
            end
        end
        # All simulations should have x == y
        @test equal_count == N
    end

    @testset "Extended jump system conservation" begin
        # SIR model built via extend should conserve population
        @variables S(t) I(t) R(t)
        @parameters β γ

        # S -> I
        j1 = ConstantRateJump(β * S * I, [S ~ Pre(S) - 1, I ~ Pre(I) + 1])
        @named sys1 = JumpSystem([j1], t, [S, I], [β])

        # I -> R
        j2 = ConstantRateJump(γ * I, [I ~ Pre(I) - 1, R ~ Pre(R) + 1])
        @named sys2 = JumpSystem([j2], t, [I, R], [γ])

        extended = extend(sys2, sys1; name = :sir)
        completed = complete(extended)

        N_pop = 1000
        u0 = [S => N_pop - 1, I => 1, R => 0]
        ps = [β => 0.3 / N_pop, γ => 0.1]

        jprob = JumpProblem(completed, [u0; ps], (0.0, 100.0); rng)
        sol = solve(jprob, SSAStepper())

        # Check conservation at multiple time points
        for i in 1:length(sol.t)
            total = sol[S][i] + sol[I][i] + sol[R][i]
            @test total == N_pop
        end
    end

    @testset "Extended poissonian jump counts" begin
        # Poisson process should have expected count ≈ λ * T
        @variables X(t)
        @parameters λ
        @poissonians dN(λ)

        @named sys1 = System(Equation[], t, [], [λ])  # Empty ODE base
        @named sys2 = System([D(X) ~ dN], t, [X], []; poissonians = [dN])

        extended = extend(sys2, sys1; name = :extended)
        compiled_sys = mtkcompile(extended)

        λ_val = 5.0
        T = 10.0
        expected_count = λ_val * T

        N = 500
        counts = zeros(N)
        for i in 1:N
            jprob = JumpProblem(compiled_sys, [X => 0, λ => λ_val], (0.0, T); rng = StableRNG(i))
            sol = solve(jprob, SSAStepper())
            counts[i] = sol[X][end]
        end

        # Mean should be close to expected
        @test abs(mean(counts) - expected_count) / expected_count < 0.1
    end
end
