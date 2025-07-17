using ModelingToolkit
using LinearSolve
using SciMLBase
using StaticArrays
using SparseArrays
using Test
using ModelingToolkit: t_nounits as t, D_nounits as D, SystemCompatibilityError

@testset "Rejects non-affine systems" begin
    @variables x y
    @mtkbuild sys = System([0 ~ x^2 + y, 0 ~ x - y])
    @test_throws SystemCompatibilityError LinearProblem(sys, nothing)
end

@variables x[1:3] [irreducible = true]
@parameters p[1:3, 1:3] q[1:3]

@mtkbuild sys = System([p * x ~ q])
# sanity check
@test length(unknowns(sys)) == length(equations(sys)) == 3
A = Float64[1 2 3; 4 3.5 1.7; 5.2 1.8 9.7]
b = Float64[2, 5, 8]
ps = [p => A, q => b]

@testset "Basics" begin
    # Ensure it works without providing `u0`
    prob = LinearProblem(sys, ps)
    @test prob.u0 === nothing
    @test SciMLBase.isinplace(prob)
    @test prob.A ≈ A
    @test prob.b ≈ b
    @test eltype(prob.A) == Float64
    @test eltype(prob.b) == Float64

    @test prob.ps[p * q] ≈ A * b

    sol = solve(prob)
    @test SciMLBase.successful_retcode(sol)
    @test prob.A * sol.u - prob.b≈zeros(3) atol=1e-10

    A2 = rand(3, 3)
    b2 = rand(3)
    @testset "remake" begin
        prob2 = remake(prob; p = [p => A2, q => b2])
        @test prob2.ps[p] ≈ A2
        @test prob2.ps[q] ≈ b2
        @test prob2.A ≈ A2
        @test prob2.b ≈ b2
    end

    prob.ps[p] = A2
    @test prob.A ≈ A2
    prob.ps[q] = b2
    @test prob.b ≈ b2
    A2[1, 1] = prob.ps[p[1, 1]] = 1.5
    @test prob.A ≈ A2
    b2[1] = prob.ps[q[1]] = 2.5
    @test prob.b ≈ b2

    @testset "expression = Val{true}" begin
        prob3e = LinearProblem(sys, ps; expression = Val{true})
        @test prob3e isa Expr
        prob3 = eval(prob3e)

        @test prob3.u0 === nothing
        @test SciMLBase.isinplace(prob3)
        @test prob3.A ≈ A
        @test prob3.b ≈ b
        @test eltype(prob3.A) == Float64
        @test eltype(prob3.b) == Float64

        @test prob3.ps[p * q] ≈ A * b

        sol = solve(prob3)
        @test SciMLBase.successful_retcode(sol)
        @test prob3.A * sol.u - prob3.b≈zeros(3) atol=1e-10
    end
end

@testset "With `u0`" begin
    prob = LinearProblem(sys, [x => ones(3); ps])
    @test prob.u0 ≈ ones(3)
    @test SciMLBase.isinplace(prob)
    @test eltype(prob.u0) == Float64

    # Observed should work
    @test prob[x[1] + x[2]] ≈ 2.0

    @testset "expression = Val{true}" begin
        prob3e = LinearProblem(sys, [x => ones(3); ps]; expression = Val{true})
        @test prob3e isa Expr
        prob3 = eval(prob3e)
        @test prob3.u0 ≈ ones(3)
        @test eltype(prob3.u0) == Float64
    end
end

@testset "SArray OOP form" begin
    prob = LinearProblem(sys, SVector{2}(ps))
    @test prob.A isa SMatrix{3, 3, Float64}
    @test prob.b isa SVector{3, Float64}
    @test !SciMLBase.isinplace(prob)
    @test prob.ps[p * q] ≈ A * b

    sol = solve(prob)
    # https://github.com/SciML/LinearSolve.jl/issues/532
    @test SciMLBase.successful_retcode(sol)
    @test prob.A * sol.u - prob.b≈zeros(3) atol=1e-10

    A2 = rand(3, 3)
    b2 = rand(3)
    @testset "remake" begin
        prob2 = remake(prob; p = [p => A2, q => b2])
        # Despite passing `Array` to `remake`
        @test prob2.A isa SMatrix{3, 3, Float64}
        @test prob2.b isa SVector{3, Float64}
        @test prob2.ps[p] ≈ A2
        @test prob2.ps[q] ≈ b2
        @test prob2.A ≈ A2
        @test prob2.b ≈ b2
    end

    @testset "expression = Val{true}" begin
        prob3e = LinearProblem(sys, SVector{2}(ps); expression = Val{true})
        @test prob3e isa Expr
        prob3 = eval(prob3e)
        @test prob3.A isa SMatrix{3, 3, Float64}
        @test prob3.b isa SVector{3, Float64}
        @test !SciMLBase.isinplace(prob3)
        @test prob3.ps[p * q] ≈ A * b

        sol = solve(prob3)
        # https://github.com/SciML/LinearSolve.jl/issues/532
        @test SciMLBase.successful_retcode(sol)
        @test prob3.A * sol.u - prob3.b≈zeros(3) atol=1e-10
    end
end

@testset "u0_constructor" begin
    prob = LinearProblem{false}(sys, ps; u0_constructor = x -> SArray{Tuple{size(x)...}}(x))
    @test prob.A isa SMatrix{3, 3, Float64}
    @test prob.b isa SVector{3, Float64}
    @test prob.ps[p * q] ≈ A * b
end

@testset "sparse form" begin
    prob = LinearProblem(sys, ps; sparse = true)
    @test issparse(prob.A)
    @test !issparse(prob.b)

    sol = solve(prob)
    # This might end up failing because of
    # https://github.com/SciML/LinearSolve.jl/issues/532
    @test SciMLBase.successful_retcode(sol)

    A2 = rand(3, 3)
    prob.ps[p] = A2
    @test prob.A ≈ A2
    b2 = rand(3)
    prob.ps[q] = b2
    @test prob.b ≈ b2

    A2 = rand(3, 3)
    b2 = rand(3)
    @testset "remake" begin
        prob2 = remake(prob; p = [p => A2, q => b2])
        @test issparse(prob2.A)
        @test !issparse(prob2.b)
        @test prob2.ps[p] ≈ A2
        @test prob2.ps[q] ≈ b2
        @test prob2.A ≈ A2
        @test prob2.b ≈ b2
    end

    @testset "expression = Val{true}" begin
        prob3e = LinearProblem(sys, ps; sparse = true, expression = Val{true})
        @test prob3e isa Expr
        prob3 = eval(prob3e)
        @test issparse(prob3.A)
        @test !issparse(prob3.b)

        sol = solve(prob3)
        # This might end up failing because of
        # https://github.com/SciML/LinearSolve.jl/issues/532
        @test SciMLBase.successful_retcode(sol)
    end
end
