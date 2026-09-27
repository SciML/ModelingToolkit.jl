using Test, ModelingToolkit, SciMLBase, OrdinaryDiffEq
using OrdinaryDiffEqBDF

@testset "Derivative initialization after tearing" begin
    @independent_variables t
    D = Differential(t)
    @variables y(t)[1:3] z(t) x(t)

    @testset "Eliminated derivatives, $specialize, guesses = $guess_derivatives" for
        specialize in (SciMLBase.AutoSpecialize, SciMLBase.FullSpecialize),
            guess_derivatives in (false, true)
        sys = complete(
            System(
                [zeros(3) ~ D(y[1:3]) + z * y[1:3], 0 ~ z - sum(y)],
                t, [y, z], []; name = :residual
            )
        )
        guesses = guess_derivatives ? [z => 0.0, D(y) => zeros(3)] : [z => 0.0]
        prob = DAEProblem{true, specialize}(
            sys, [y => [1.0, 2.0, 3.0], D(z) => 0.0], (0.0, 1.0); guesses
        )
        @test isempty(unknowns(prob.f.initializeprob.f.sys))
        for initkwargs in ((;), (; initializealg = SciMLBase.OverrideInit()))
            integ = init(prob, DFBDF(); initkwargs...)
            @test integ.u ≈ [1.0, 2.0, 3.0, 6.0]
            @test integ.du ≈ [-6.0, -12.0, -18.0, 0.0]
            residual = similar(integ.u)
            prob.f(residual, integ.du, integ.u, integ.p, 0.0)
            @test residual ≈ zeros(4)
        end
    end

    @testset "Derivative root, sign = $sign, system guess = $system_guess" for
        sign in (-1.0, 1.0), system_guess in (false, true)
        derivative_guesses = [D(x) => sign]
        sys = mtkcompile(
            System(
                [D(D(x)) ~ 0], t;
                initialization_eqs = [D(x)^2 ~ 1, x ~ 0],
                guesses = system_guess ? derivative_guesses : [], name = :second_order
            )
        )
        prob = ODEProblem(
            sys, [], (0.0, 1.0); guesses = system_guess ? [] : derivative_guesses
        )
        @test prob.f.initializeprob.u0 ≈ [sign]
        sol = solve(prob, Tsit5())
        @test SciMLBase.successful_retcode(sol)
        @test sol(1.0; idxs = x) ≈ sign atol = 1.0e-6
    end
end
