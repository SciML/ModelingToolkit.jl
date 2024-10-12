using ModelingToolkit, NonlinearSolve, SymbolicIndexingInterface
using LinearAlgebra
using Test
import HomotopyContinuation

@testset "No parameters" begin
    @variables x y z
    eqs = [0 ~ x^2 + y^2 + 2x * y
           0 ~ x^2 + 4x + 4
           0 ~ y * z + 4x^2]
    @mtkbuild sys = NonlinearSystem(eqs)
    prob = HomotopyContinuationProblem(sys, [x => 1.0, y => 1.0, z => 1.0], [])
    @test prob[x] == prob[y] == prob[z] == 1.0
    @test prob[x + y] == 2.0
    sol = solve(prob; threading = false)
    @test SciMLBase.successful_retcode(sol)
    @test norm(sol.resid)≈0.0 atol=1e-10
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

    @mtkbuild sys = NonlinearSystem(eqs)
    prob = HomotopyContinuationProblem(sys, [x => 1.0, y => 1.0, z => 1.0],
        [p => 2.0, q => 4, r => Wrapper([1.0 1.0; 0.0 0.0])])
    @test prob.ps[p] == 2.0
    @test prob.ps[q] == 4
    @test prob.ps[r].x == [1.0 1.0; 0.0 0.0]
    @test prob.ps[p * q] == 8.0
    sol = solve(prob; threading = false)
    @test SciMLBase.successful_retcode(sol)
    @test norm(sol.resid)≈0.0 atol=1e-10
end

@testset "Array variables" begin
    @variables x[1:3]
    @parameters p[1:3]
    _x = collect(x)
    eqs = collect(0 .~ vec(sum(_x * _x'; dims = 2)) + collect(p))
    @mtkbuild sys = NonlinearSystem(eqs)
    prob = HomotopyContinuationProblem(sys, [x => ones(3)], [p => 1:3])
    @test prob[x] == ones(3)
    @test prob[p + x] == [2, 3, 4]
    prob[x] = 2ones(3)
    @test prob[x] == 2ones(3)
    prob.ps[p] = [2, 3, 4]
    @test prob.ps[p] == [2, 3, 4]
    sol = @test_nowarn solve(prob; threading = false)
    @test sol.retcode == ReturnCode.ConvergenceFailure
end
