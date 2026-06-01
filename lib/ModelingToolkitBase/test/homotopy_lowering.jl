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

    @testset "L1-Q5 OverrideInitData metadata exposes setp handle + default alg" begin
        using ModelingToolkitBase: System, mtkcompile, HomotopySweep, TrivialThenSweep
        using ModelingToolkitBase: t_nounits as t, D_nounits as D
        using OrdinaryDiffEqRosenbrock: Rodas5P
        using OrdinaryDiffEqNonlinearSolve
        using SymbolicIndexingInterface
        using SciMLBase

        @variables x(t) y(t)
        @parameters p
        eqs = [D(x) ~ -x,
               0 ~ homotopy(y^2 - p, y - 1)]
        @named sys = System(eqs, t; guesses = [y => 2.5])
        sys = mtkcompile(sys)
        prob = ODEProblem(sys, Dict(x => 1.0, p => 9.0), (0.0, 1.0))

        meta = prob.f.initialization_data.metadata
        @test meta.homotopy_set_λ! !== nothing
        @test meta.homotopy_default_initializealg !== nothing

        # setp handle round-trip — write λ = 0.3, read it back
        initprob = prob.f.initialization_data.initializeprob
        set_λ! = meta.homotopy_set_λ!
        new_p = set_λ!(parameter_values(initprob), 0.3)
        getter = SymbolicIndexingInterface.getp(initprob, :__homotopy_λ)
        @test isapprox(getter(new_p), 0.3; atol = 1e-12)

        # Default initializealg shape
        default_alg = meta.homotopy_default_initializealg
        @test default_alg isa SciMLBase.OverrideInit
        @test default_alg.nlsolve isa TrivialThenSweep
        @test default_alg.nlsolve.sweep isa HomotopySweep

        # Non-homotopy system: metadata fields stay `nothing`.
        @variables a(t)
        @parameters q
        @named sys2 = System([D(a) ~ -q * a], t; guesses = [a => 1.0])
        sys2 = mtkcompile(sys2)
        prob2 = ODEProblem(sys2, Dict(a => 1.0, q => 2.0), (0.0, 1.0))
        if prob2.f.initialization_data !== nothing
            meta2 = prob2.f.initialization_data.metadata
            @test meta2.homotopy_set_λ! === nothing
            @test meta2.homotopy_default_initializealg === nothing
        end
    end

    @testset "L1-Q6 default initializealg injected into prob.kwargs for homotopy systems" begin
        using ModelingToolkitBase: System, mtkcompile, TrivialThenSweep, TrivialHomotopy
        using ModelingToolkitBase: t_nounits as t, D_nounits as D
        using OrdinaryDiffEqRosenbrock
        using OrdinaryDiffEqNonlinearSolve
        using SciMLBase

        @variables x(t) y(t)
        @parameters p
        eqs = [D(x) ~ -x,
               0 ~ homotopy(y^2 - p, y - 1)]
        @named sys = System(eqs, t; guesses = [y => 2.5])
        sys = mtkcompile(sys)

        # Auto-injection path: user passes no `initializealg`
        prob_default = ODEProblem(sys, Dict(x => 1.0, p => 9.0), (0.0, 1.0))
        @test haskey(prob_default.kwargs, :initializealg)
        injected = prob_default.kwargs[:initializealg]
        @test injected isa SciMLBase.OverrideInit
        @test injected.nlsolve isa TrivialThenSweep

        # Explicit override path: user passes their own initializealg
        explicit_alg = SciMLBase.OverrideInit(; nlsolve = TrivialHomotopy())
        prob_explicit = ODEProblem(sys, Dict(x => 1.0, p => 9.0), (0.0, 1.0);
                                    initializealg = explicit_alg)
        @test haskey(prob_explicit.kwargs, :initializealg)
        # User's explicit alg must win over MTK's default
        @test prob_explicit.kwargs[:initializealg] === explicit_alg

        # Non-homotopy system: no `initializealg` injection
        @variables a(t)
        @parameters q
        @named sys2 = System([D(a) ~ -q * a], t; guesses = [a => 1.0])
        sys2 = mtkcompile(sys2)
        prob_nohomotopy = ODEProblem(sys2, Dict(a => 1.0, q => 2.0), (0.0, 1.0))
        @test !haskey(prob_nohomotopy.kwargs, :initializealg)
    end

    @testset "L1-Q7 observed equations are homotopy-free after lowering" begin
        # Regression guard: `add_homotopy_parameter` lowers homotopy nodes in
        # observed equations too (PressureDrop-style: an eliminated variable's
        # definition lives in observed). Downstream observed codegen must never
        # see an opaque `homotopy(...)` node — assert observed is homotopy-free
        # after lowering a System whose homotopy-defined variable is eliminated.
        using ModelingToolkitBase: System, mtkcompile, observed, has_homotopy
        using ModelingToolkitBase: t_nounits as t, D_nounits as D
        using OrdinaryDiffEqRosenbrock
        using OrdinaryDiffEqNonlinearSolve

        # `w ~ homotopy(...)` defines w explicitly from a state, so mtkcompile
        # eliminates w into observed — the codepath where an eliminated
        # variable's definition (PressureDrop's m_flow) lands in observed.
        @variables x(t) w(t)
        @parameters p
        eqs = [D(x) ~ -x,
               w ~ homotopy(x^2 - p, x - 1)]
        @named sys = System(eqs, t; guesses = [w => 2.5])
        sys = mtkcompile(sys)

        obs = observed(sys)
        @test obs !== nothing && !isempty(obs)
        for eq in obs
            @test !has_homotopy(eq.lhs)
            @test !has_homotopy(eq.rhs)
        end

        # The lowered ODEProblem's observed must likewise be homotopy-free.
        prob = ODEProblem(sys, Dict(x => 1.0, p => 9.0), (0.0, 1.0))
        for eq in observed(prob.f.sys)
            @test !has_homotopy(eq.lhs)
            @test !has_homotopy(eq.rhs)
        end
    end
end
