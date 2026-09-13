using ModelingToolkitBase, Symbolics, SciMLBase, OrdinaryDiffEqBDF, SparseArrays, LinearAlgebra, Test
using DiffEqBase: BrownFullBasicInit

@testset "Implicit array residual Jacobians include both sides" begin
    @independent_variables t
    @variables x(t)[1:3]
    D = Differential(t)
    sys = complete(System([D(x[1:2]) ~ -x[1:2], x[3] ~ x[1] + 2x[2]], t, collect(x), []; name = :array_dae))
    op = [x => [1.0, 2.0, 5.0], D(x) => [-1.0, -2.0, 0.0]]
    for sparse in (false, true)
        prob = DAEProblem(sys, op, (0.0, 1.0); jac = true, sparse, build_initializeprob = false)
        matrix = sparse ? copy(prob.f.jac_prototype) : zeros(3, 3)
        prob.f.jac(matrix, prob.du0, prob.u0, prob.p, 4.0, 0.0)
        @test Matrix(matrix) ≈ [-5 0 0;0 -5 0;1 2 -1]
        sol = solve(prob, DFBDF(); reltol = 1.0e-8, abstol = 1.0e-10, initializealg = BrownFullBasicInit())
        @test SciMLBase.successful_retcode(sol)
        @test sol.u[end] ≈ exp(-1) .* [1.0, 2.0, 5.0] rtol = 1.0e-6
    end
    prototype = DAEFunction(sys; sparse = true).jac_prototype
    @test size(prototype) == (3, 3)
    @test nnz(prototype) == 5
end

@testset "DAE prototype retains cancelling state and derivative entries" begin
    @independent_variables t
    @variables x(t)
    D = Differential(t)
    sys = complete(System([D(x) ~ x], t, [x], []; name = :decay))
    f = DAEFunction(sys; jac = true, sparse = true)
    @test nnz(f.jac_prototype) == 1
    J = copy(f.jac_prototype)
    for gamma in (1.0, 2.0)
        f.jac(J, [1.0], [1.0], SciMLBase.NullParameters(), gamma, 0.0)
        @test J[1, 1] == 1 - gamma
    end
end
@testset "Empty array residuals add no Jacobian rows" begin
    @independent_variables t
    @variables x(t)[1:3]
    D = Differential(t)
    empty_lhs = Symbolics.term(
        +, Symbolics.unwrap(x[2:1]), Symbolics.unwrap(x[3:2]);
        type = Vector{Float64}, shape = (1:0,)
    )
    sys = complete(System([D(x) ~ -x, empty_lhs ~ zeros(0)], t, collect(x), []; name = :empty_residual))
    f = DAEFunction(sys; sparse = true, jac = true)
    @test size(f.jac_prototype) == (3, 3)
    J = copy(f.jac_prototype)
    f.jac(J, -ones(3), ones(3), SciMLBase.NullParameters(), 2.0, 0.0)
    @test Matrix(J) == -3 .* Matrix{Float64}(I, 3, 3)
end
