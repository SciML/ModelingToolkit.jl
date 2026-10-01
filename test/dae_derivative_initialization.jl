using Test, ModelingToolkit, SciMLBase, OrdinaryDiffEq
using OrdinaryDiffEqBDF

@testset "Derivative initialization after tearing" begin
    @independent_variables t
    D = Differential(t)
    @variables y(t)[1:3] z(t) x(t) x1(t) x2(t) w(t)
    @parameters k
    array_sys = complete(
        System(
            [zeros(3) ~ D(y[1:3]) + z * y[1:3], 0 ~ z - sum(y)], t, [y, z], [];
            name = :residual
        )
    )
    scalar_sys = complete(
        System(
            [0 ~ D(x1) + k * w * x1, 0 ~ D(x2) + k * w * x2, 0 ~ w - x1 - x2],
            t, [x1, x2, w], [k]; name = :scalar
        )
    )

    @testset "Eliminated derivatives, $specialize, guesses = $guess_derivatives" for
        specialize in (SciMLBase.AutoSpecialize, SciMLBase.FullSpecialize),
            guess_derivatives in (false, true)
        guesses = guess_derivatives ? [z => 0.0, D(y) => zeros(3)] : [z => 0.0]
        prob = DAEProblem{true, specialize}(
            array_sys, [y => [1.0, 2.0, 3.0], D(z) => 0.0], (0.0, 1.0); guesses
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

    @testset "Omitted derivatives give a consistent du0: $label" for (label, sys, op, guesses, du0) in (
            (
                "array, D(z) given", array_sys,
                [y => [1.0, 2.0, 3.0], D(z) => 0.0], [z => 0.0], [-6.0, -12.0, -18.0, 0.0],
            ),
            ("array, none given", array_sys, [y => [1.0, 2.0, 3.0]], [z => 0.0], [-6.0, -12.0, -18.0, 0.0]),
            ("scalar, none given", scalar_sys, [x1 => 1.0, x2 => 2.0, k => 1.0], [w => 0.0], [-3.0, -6.0, 0.0]),
            (
                "scalar, D(w) given", scalar_sys,
                [x1 => 1.0, x2 => 2.0, k => 1.0, D(w) => 0.0], [w => 0.0], [-3.0, -6.0, 0.0],
            ),
        )
        prob = DAEProblem(sys, op, (0.0, 1.0); guesses)
        @test prob.du0 ≈ du0
        residual = similar(prob.u0)
        prob.f(residual, prob.du0, prob.u0, prob.p, 0.0)
        @test residual ≈ zeros(length(residual)) atol = 1.0e-12
        sol = solve(prob, DFBDF(); initializealg = SciMLBase.CheckInit())
        @test SciMLBase.successful_retcode(sol)
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
