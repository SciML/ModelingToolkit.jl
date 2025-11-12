using ModelingToolkit
using OrdinaryDiffEq
using LinearAlgebra
using Test
using ModelingToolkit: t_nounits as t, D_nounits as D

# from https://docs.sciml.ai/SciMLBenchmarksOutput/dev/AstroChem/nelson/
@testset "Astrochem model" begin
    function Nelson!(du, u, p, t)
        T, Av, Go, n_H, shield = p

        # 1: H2
        du[1] = -1.2e-17 * u[1] +
                n_H * (1.9e-6 * u[2] * u[3]) / (T^0.54) -
                n_H * 4e-16 * u[1] * u[12] -
                n_H * 7e-15 * u[1] * u[5] +
                n_H * 1.7e-9 * u[10] * u[2] +
                n_H * 2e-9 * u[2] * u[6] +
                n_H * 2e-9 * u[2] * u[14] +
                n_H * 8e-10 * u[2] * u[8] + sin(u[1]) / 1e20 # artificial nonlinear term for testing

        # 2: H3+
        du[2] = 1.2e-17 * u[1] +
                n_H * (-1.9e-6 * u[3] * u[2]) / (T^0.54) -
                n_H * 1.7e-9 * u[10] * u[2] -
                n_H * 2e-9 * u[2] * u[6] -
                n_H * 2e-9 * u[2] * u[14] -
                n_H * 8e-10 * u[2] * u[8]

        # 3: e
        du[3] = n_H * (-1.4e-10 * u[3] * u[12]) / (T^0.61) -
                n_H * (3.8e-10 * u[13] * u[3]) / (T^0.65) -
                n_H * (3.3e-5 * u[11] * u[3]) / T +
                1.2e-17 * u[1] -
                n_H * (1.9e-6 * u[3] * u[2]) / (T^0.54) +
                6.8e-18 * u[4] -
                n_H * (9e-11 * u[3] * u[5]) / (T^0.64) +
                3e-10 * Go * exp(-3 * Av) * u[6] +
                n_H * 2e-9 * u[2] * u[13]
        +2.0e-10 * Go * exp(-1.9 * Av) * u[14]

        # 4: He
        du[4] = n_H * (9e-11 * u[3] * u[5]) / (T^0.64) -
                6.8e-18 * u[4] +
                n_H * 7e-15 * u[1] * u[5] +
                n_H * 1.6e-9 * u[10] * u[5]

        # 5: He+
        du[5] = 6.8e-18 * u[4] -
                n_H * (9e-11 * u[3] * u[5]) / (T^0.64) -
                n_H * 7e-15 * u[1] * u[5] -
                n_H * 1.6e-9 * u[10] * u[5]

        # 6: C
        du[6] = n_H * (1.4e-10 * u[3] * u[12]) / (T^0.61) -
                n_H * 2e-9 * u[2] * u[6] -
                n_H * 5.8e-12 * (T^0.5) * u[9] * u[6] +
                1e-9 * Go * exp(-1.5 * Av) * u[7] -
                3e-10 * Go * exp(-3 * Av) * u[6] +
                1e-10 * Go * exp(-3 * Av) * u[10] * shield

        # 7: CHx
        du[7] = n_H * (-2e-10) * u[7] * u[8] +
                n_H * 4e-16 * u[1] * u[12] +
                n_H * 2e-9 * u[2] * u[6] -
                1e-9 * Go * u[7] * exp(-1.5 * Av)

        # 8: O
        du[8] = n_H * (-2e-10) * u[7] * u[8] +
                n_H * 1.6e-9 * u[10] * u[5] -
                n_H * 8e-10 * u[2] * u[8] +
                5e-10 * Go * exp(-1.7 * Av) * u[9] +
                1e-10 * Go * exp(-3 * Av) * u[10] * shield

        # 9: OHx
        du[9] = n_H * (-1e-9) * u[9] * u[12] +
                n_H * 8e-10 * u[2] * u[8] -
                n_H * 5.8e-12 * (T^0.5) * u[9] * u[6] -
                5e-10 * Go * exp(-1.7 * Av) * u[9]

        # 10: CO
        du[10] = n_H * (3.3e-5 * u[11] * u[3]) / T +
                 n_H * 2e-10 * u[7] * u[8] -
                 n_H * 1.7e-9 * u[10] * u[2] -
                 n_H * 1.6e-9 * u[10] * u[5] +
                 n_H * 5.8e-12 * (T^0.5) * u[9] * u[6] -
                 1e-10 * Go * exp(-3 * Av) * u[10] +
                 1.5e-10 * Go * exp(-2.5 * Av) * u[11] * shield

        # 11: HCO+
        du[11] = n_H * (-3.3e-5 * u[11] * u[3]) / T +
                 n_H * 1e-9 * u[9] * u[12] +
                 n_H * 1.7e-9 * u[10] * u[2] -
                 1.5e-10 * Go * exp(-2.5 * Av) * u[11]

        # 12: C+
        du[12] = n_H * (-1.4e-10 * u[3] * u[12]) / (T^0.61) -
                 n_H * 4e-16 * u[1] * u[12] -
                 n_H * 1e-9 * u[9] * u[12] +
                 n_H * 1.6e-9 * u[10] * u[5] +
                 3e-10 * Go * exp(-3 * Av) * u[6]

        # 13: M+
        du[13] = n_H * (-3.8e-10 * u[13] * u[3]) / (T^0.65) +
                 n_H * 2e-9 * u[2] * u[14] +
                 2.0e-10 * Go * exp(-1.9 * Av) * u[14]

        # 14: M
        du[14] = n_H * (3.8e-10 * u[13] * u[3]) / (T^0.65) -
                 n_H * 2e-9 * u[2] * u[14] -
                 2.0e-10 * Go * exp(-1.9 * Av) * u[14]
    end

    # Set the Timespan, Parameters, and Initial Conditions
    seconds_per_year = 3600 * 24 * 365
    tspan = (0.0, 30000 * seconds_per_year) # ~30 thousand yrs

    params = (10,  # T
        2,   # Av
        1.7, # Go
        611, # n_H
        1)   # shield

    u0 = [0.5,      # 1:  H2
        9.059e-9, # 2:  H3+
        2.0e-4,   # 3:  e
        0.1,      # 4:  He
        7.866e-7, # 5:  He+
        0.0,      # 6:  C
        0.0,      # 7:  CHx
        0.0004,   # 8:  O
        0.0,      # 9:  OHx
        0.0,      # 10: CO
        0.0,      # 11: HCO+
        0.0002,   # 12: C+
        2.0e-7,   # 13: M+
        2.0e-7]   # 14: M

    prob = ODEProblem(Nelson!, u0, tspan, params)
    sys = mtkcompile(modelingtoolkitize(prob))
    A, B, C = ModelingToolkit.calculate_semiquadratic_form(sys)
    @test A !== nothing
    @test B !== nothing
    @test C !== nothing
    x = unknowns(sys)
    linear_expr = A * x
    linear_fun, = generate_custom_function(sys, linear_expr; expression = Val{false})
    quadratic_expr = reduce(vcat, [x' * _B * x for _B in B if _B !== nothing])
    quadratic_fun, = generate_custom_function(sys, quadratic_expr; expression = Val{false})
    nonlinear_expr = C
    nonlinear_fun, = generate_custom_function(sys, nonlinear_expr; expression = Val{false})
    prob = ODEProblem(sys, nothing, tspan)
    linear_val = linear_fun(prob.u0, prob.p, 0.0)
    quadratic_val = quadratic_fun(prob.u0, prob.p, 0.0)
    nonlinear_val = nonlinear_fun(prob.u0, prob.p, 0.0)
    refsol = solve(prob, Vern9(); abstol = 1e-14, reltol = 1e-14)

    @testset "stiff_linear: $stiff_linear, stiff_quadratic: $stiff_quadratic, stiff_nonlinear: $stiff_nonlinear" for (
        stiff_linear, stiff_quadratic, stiff_nonlinear) in Iterators.product(
        [false, true], [false, true], [false, true])
        kwargs = (; stiff_linear, stiff_quadratic, stiff_nonlinear)
        if stiff_linear == stiff_quadratic == stiff_nonlinear
            if stiff_linear
                @test_throws ["All of", "cannot be stiff"] SemilinearODEProblem(
                    sys, nothing, tspan; kwargs...)
                @test_throws ["All of", "cannot be stiff"] SemilinearODEFunction(
                    sys; kwargs...)
            else
                @test_throws ["All of", "cannot be non-stiff"] SemilinearODEProblem(
                    sys, nothing, tspan; kwargs...)
                @test_throws ["All of", "cannot be non-stiff"] SemilinearODEFunction(
                    sys; kwargs...)
            end
            continue
        end

        reference_f1 = zeros(length(u0))
        reference_f2 = zeros(length(u0))
        mul!(stiff_linear ? reference_f1 : reference_f2, I, linear_val, true, true)
        mul!(stiff_quadratic ? reference_f1 : reference_f2, I, quadratic_val, true, true)
        mul!(stiff_nonlinear ? reference_f1 : reference_f2, I, nonlinear_val, true, true)

        @testset "Standard" begin
            prob = SemilinearODEProblem(sys, nothing, tspan; kwargs...)
            @test prob.f.f1(prob.u0, prob.p, 0.0)≈reference_f1 atol=1e-10
            @test prob.f.f2(prob.u0, prob.p, 0.0)≈reference_f2 atol=1e-10
            sol = solve(prob, KenCarp47())
            @test SciMLBase.successful_retcode(sol)
            @test refsol(sol.t).u≈sol.u atol=1e-8 rtol=1e-8
        end

        @testset "Symbolic jacobian" begin
            prob = SemilinearODEProblem(sys, nothing, tspan; jac = true, kwargs...)
            @test prob.f.f1.jac !== nothing
            sol = solve(prob, KenCarp47())
            @test SciMLBase.successful_retcode(sol)
            @test refsol(sol.t).u≈sol.u atol=1e-8 rtol=1e-8
        end

        @testset "Sparse" begin
            @test_throws ["sparse form", "unavailable"] SemilinearODEProblem(
                sys, nothing, tspan; sparse = true, kwargs...)
            @test_skip begin
                sol = solve(prob, KenCarp47())
                @test SciMLBase.successful_retcode(sol)
                @test refsol(sol.t).u≈sol.u atol=1e-8 rtol=1e-8
            end
        end

        @testset "Sparsejac" begin
            @test_throws ["sparse form", "unavailable"] SemilinearODEProblem(
                sys, nothing, tspan; jac = true, sparse = true, kwargs...)
        end
    end
end
