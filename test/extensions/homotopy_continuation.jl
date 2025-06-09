using ModelingToolkit, NonlinearSolve, NonlinearSolveHomotopyContinuation,
      SymbolicIndexingInterface
using SymbolicUtils
import ModelingToolkit as MTK
using LinearAlgebra
using Test

allrootsalg = HomotopyContinuationJL{true}(; threading = false)
singlerootalg = HomotopyContinuationJL{false}(; threading = false)

function test_single_root(sol; atol = 1e-10)
    @test SciMLBase.successful_retcode(sol)
    @test norm(sol.resid)≈0.0 atol=atol
end

function test_all_roots(sol; atol = 1e-10)
    @test sol.converged
    for nlsol in sol.u
        @test SciMLBase.successful_retcode(nlsol)
        @test norm(nlsol.resid)≈0.0 atol=1e-10
    end
end

function solve_allroots_closest(prob)
    sol = solve(prob, allrootsalg)
    return argmin(sol.u) do nlsol
        return norm(nlsol.u - prob.u0)
    end
end

@testset "No parameters" begin
    @variables x y z
    eqs = [0 ~ x^2 + y^2 + 2x * y
           0 ~ x^2 + 4x + 4
           0 ~ y * z + 4x^2]
    @mtkcompile sys = System(eqs)
    u0 = [x => 1.0, y => 1.0, z => 1.0]
    prob = HomotopyContinuationProblem(sys, u0)
    @test prob isa NonlinearProblem
    @test prob[x] == prob[y] == prob[z] == 1.0
    @test prob[x + y] == 2.0
    sol = solve(prob, singlerootalg)
    test_single_root(sol)
    sol = solve(prob, allrootsalg)
    test_all_roots(sol)
end

struct Wrapper
    x::Matrix{Float64}
end

@testset "Parameters" begin
    wrapper(w::Wrapper) = det(w.x)
    @register_symbolic wrapper(w::Wrapper)

    @variables x y z
    @parameters p q::Int r::Wrapper

    eqs = [0 ~ x^2 + y^2 + p * x * y
           0 ~ x^2 + 4x + q
           0 ~ y * z + 4x^2 + wrapper(r)]

    @mtkcompile sys = System(eqs)
    prob = HomotopyContinuationProblem(sys,
        [x => 1.0, y => 1.0, z => 1.0, p => 2.0, q => 4, r => Wrapper([1.0 1.0; 0.0 0.0])])
    @test prob.ps[p] == 2.0
    @test prob.ps[q] == 4
    @test prob.ps[r].x == [1.0 1.0; 0.0 0.0]
    @test prob.ps[p * q] == 8.0
    sol = solve(prob, singlerootalg)
    test_single_root(sol)
    sol = solve(prob, allrootsalg)
    test_all_roots(sol)
end

@testset "Array variables" begin
    @variables x[1:3]
    @parameters p[1:3]
    _x = collect(x)
    eqs = collect(0 .~ vec(sum(_x * _x'; dims = 2)) + collect(p))
    @mtkcompile sys = System(eqs)
    prob = HomotopyContinuationProblem(sys, [x => ones(3), p => 1:3])
    @test prob[x] == ones(3)
    @test prob[p + x] == [2, 3, 4]
    prob[x] = 2ones(3)
    @test prob[x] == 2ones(3)
    prob.ps[p] = [2, 3, 4]
    @test prob.ps[p] == [2, 3, 4]
    sol = @test_nowarn solve(prob, singlerootalg)
    @test sol.retcode == ReturnCode.ConvergenceFailure
end

@testset "Parametric exponents" begin
    @variables x = 1.0
    @parameters n::Integer = 4
    @mtkcompile sys = System([x^n + x^2 - 1 ~ 0])
    prob = HomotopyContinuationProblem(sys, [])
    sol = solve(prob, singlerootalg)
    test_single_root(sol)
    sol = solve(prob, allrootsalg)
    test_all_roots(sol)
end

@testset "Polynomial check and warnings" begin
    @variables x = 1.0
    @mtkcompile sys = System([x^1.5 + x^2 - 1 ~ 0])
    @test_throws ["Cannot convert", "Unable", "symbolically solve",
        "Exponent", "not an integer", "not a polynomial"] HomotopyContinuationProblem(
        sys, [])

    @mtkcompile sys = System([x^x - x ~ 0])
    @test_throws ["Cannot convert", "Unable", "symbolically solve",
        "Exponent", "unknowns", "not a polynomial"] HomotopyContinuationProblem(
        sys, [])
    @mtkcompile sys = System([((x^2) / sin(x))^2 + x ~ 0])
    @test_throws ["Cannot convert", "both polynomial", "non-polynomial",
        "recognized", "sin", "not a polynomial"] HomotopyContinuationProblem(
        sys, [])

    @variables y = 2.0
    @mtkcompile sys = System([x^2 + y^2 + 2 ~ 0, y ~ sin(x)])
    @test_throws ["Cannot convert", "recognized", "sin", "not a polynomial"] HomotopyContinuationProblem(
        sys, [])

    @mtkcompile sys = System([x^2 + y^2 - 2 ~ 0, sin(x + y) ~ 0])
    @test_throws ["Cannot convert", "function of multiple unknowns"] HomotopyContinuationProblem(
        sys, [])

    @mtkcompile sys = System([sin(x)^2 + 1 ~ 0, cos(y) - cos(x) - 1 ~ 0])
    @test_throws ["Cannot convert", "multiple non-polynomial terms", "same unknown"] HomotopyContinuationProblem(
        sys, [])

    @mtkcompile sys = System([sin(x^2)^2 + sin(x^2) - 1 ~ 0])
    @test_throws ["import Nemo"] HomotopyContinuationProblem(sys, [])
end

import Nemo

@testset "With Nemo" begin
    @variables x = 2.0
    @mtkcompile sys = System([sin(x^2)^2 + sin(x^2) - 1 ~ 0])
    prob = HomotopyContinuationProblem(sys, [])
    @test prob[1] ≈ 2.0
    # singlerootalg doesn't converge
    sol = solve(prob, allrootsalg).u[1]
    _x = sol[1]
    @test sin(_x^2)^2 + sin(_x^2) - 1≈0.0 atol=1e-12
end

@testset "Function of polynomial" begin
    @variables x=0.25 y=0.125
    a = sin(x^2 - 4x + 1)
    b = cos(3log(y) + 4)
    @mtkcompile sys = System([(a^2 - 5a * b + 6b^2) / (a - 0.25) ~ 0
                              (a^2 - 0.75a + 0.125) ~ 0])
    prob = HomotopyContinuationProblem(sys, [])
    @test prob[x] ≈ 0.25
    @test prob[y] ≈ 0.125
    sol = solve(prob, allrootsalg).u[1]
    @test SciMLBase.successful_retcode(sol)
    @test sol[a]≈0.5 atol=1e-6
    @test isapprox(sol[b], 0.25; atol = 1e-6) || isapprox(sol[b], 0.5 / 3; atol = 1e-6)
end

@testset "Rational functions" begin
    @variables x=2.0 y=2.0
    @parameters n = 5
    @mtkcompile sys = System([
        0 ~ (x^2 - n * x + 6) * (x - 1) / (x - 2) / (x - 3)
    ])
    prob = HomotopyContinuationProblem(sys, [])
    sol = solve_allroots_closest(prob)
    @test sol[x] ≈ 1.0
    p = parameter_values(prob)
    for invalid in [2.0, 3.0]
        for err in [-9e-8, 0, 9e-8]
            @test any(<=(1e-7), prob.f.denominator([invalid + err, 2.0], p))
        end
    end

    @named sys = System(
        [
            0 ~ (x - 2) / (x - 4) + ((x - 3) / (y - 7)) / ((x^2 - 4x + y) / (x - 2.5)),
            0 ~ ((y - 3) / (y - 4)) * (n / (y - 5)) + ((x - 1.5) / (x - 5.5))^2
        ],
        [x, y],
        [n]; defaults = [n => 4])
    sys = complete(sys)
    prob = HomotopyContinuationProblem(sys, [])
    sol = solve(prob, singlerootalg)
    disallowed_x = [4, 5.5]
    disallowed_y = [7, 5, 4]
    @test all(!isapprox(sol[x]; atol = 1e-8), disallowed_x)
    @test all(!isapprox(sol[y]; atol = 1e-8), disallowed_y)
    @test abs(sol[x ^ 2 - 4x + y]) >= 1e-8

    p = parameter_values(prob)
    for val in disallowed_x
        for err in [-9e-8, 0, 9e-8]
            @test any(<=(1e-7), prob.f.denominator([val + err, 2.0], p))
        end
    end
    for val in disallowed_y
        for err in [-9e-8, 0, 9e-8]
            @test any(<=(1e-7), prob.f.denominator([2.0, val + err], p))
        end
    end
    @test prob.f.denominator([2.0, 4.0], p)[1] <= 1e-8

    @testset "Rational function in observed" begin
        @variables x=1 y=1
        @mtkcompile sys = System([x^2 + y^2 - 2x - 2 ~ 0, y ~ (x - 1) / (x - 2)])
        prob = HomotopyContinuationProblem(sys, [])
        @test any(prob.f.denominator([2.0], parameter_values(prob)) .≈ 0.0)
        @test SciMLBase.successful_retcode(solve(prob, singlerootalg))
    end

    @testset "Rational function forced to common denominators" begin
        @variables x = 1
        @mtkcompile sys = System([0 ~ 1 / (1 + x) - x])
        prob = HomotopyContinuationProblem(sys, [])
        @test any(prob.f.denominator([-1.0], parameter_values(prob)) .≈ 0.0)
        sol = solve(prob, singlerootalg)
        @test SciMLBase.successful_retcode(sol)
        @test 1 / (1 + sol.u[1]) - sol.u[1]≈0.0 atol=1e-10
    end
end

@testset "Non-polynomial observed not used in equations" begin
    @variables x=1 y
    @mtkcompile sys = System([x^2 - 2 ~ 0, y ~ sin(x)])
    prob = HomotopyContinuationProblem(sys, [])
    sol = @test_nowarn solve(prob, singlerootalg)
    @test sol[x] ≈ √2.0
    @test sol[y] ≈ sin(√2.0)
end

@testset "`fraction_cancel_fn`" begin
    @variables x = 1
    @named sys = System([0 ~ ((x^2 - 5x + 6) / (x - 2) - 1) * (x^2 - 7x + 12) /
                             (x - 4)^3])
    sys = complete(sys)

    @testset "`simplify_fractions`" begin
        prob = HomotopyContinuationProblem(sys, [])
        @test prob.f.denominator([0.0], parameter_values(prob)) ≈ [4.0]
    end
    @testset "`nothing`" begin
        prob = HomotopyContinuationProblem(sys, []; fraction_cancel_fn = nothing)
        @test sort(prob.f.denominator([0.0], parameter_values(prob))) ≈ [2.0, 4.0^3]
    end
end
