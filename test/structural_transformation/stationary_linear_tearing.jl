using Test
using ModelingToolkit
using NonlinearSolve

# 5-point Laplace stencil on an n x n grid with u = x * y on the boundary, scaled by 1 / h^2
# so that the centre coefficient exceeds the integer-linear subsystem's coefficient cap.
function laplace_system(n; nonlinear = false, coupled = false)
    h = 1 / (n - 1)
    @variables u[1:n, 1:n]
    eqs = Equation[]
    for i in 1:n, j in 1:n
        x, y = (i - 1) * h, (j - 1) * h
        if i in (1, n) || j in (1, n)
            push!(eqs, u[i, j] ~ x * y)
        else
            stencil = (u[i + 1, j] + u[i - 1, j] + u[i, j + 1] + u[i, j - 1] - 4 * u[i, j]) / h^2
            push!(
                eqs, 0 ~ stencil + (nonlinear ? 0.01 * u[i, j]^2 - 0.01 * x^2 * y^2 : 0) +
                    (coupled ? 0.01 * u[i, j] * u[i + 1, j] : 0)
            )
        end
    end
    return System(eqs, nothing, vec(collect(u)), []; name = :laplace), u
end

max_error(sol, u, n) = maximum(abs(sol[u[i, j]] - (i - 1) * (j - 1) / (n - 1)^2) for i in 1:n, j in 1:n)

@testset "Large linear blocks of stationary systems are not torn" begin
    for n in (8, 9, 11)
        sys, u = laplace_system(n)
        sys = mtkcompile(sys)
        @test length(unknowns(sys)) == (n - 2)^2
        prob = NonlinearProblem(sys, Dict(vec(collect(u)) .=> 0.0))
        sol = solve(prob, NewtonRaphson())
        @test SciMLBase.successful_retcode(sol)
        @test max_error(sol, u, n) < 1.0e-12
    end

    # The previous behavior is available.
    sys, u = laplace_system(9)
    @test length(unknowns(mtkcompile(sys; max_linear_scc_tearing = typemax(Int)))) < 7^2
end

@testset "Small linear blocks and nonlinear blocks are still torn" begin
    # 16 or fewer equations in the block
    sys, _ = laplace_system(6)
    @test length(unknowns(mtkcompile(sys))) == length(unknowns(mtkcompile(sys; max_linear_scc_tearing = typemax(Int))))
    @test length(unknowns(mtkcompile(sys))) < 4^2

    # The same grid with a nonlinear term
    sys, u = laplace_system(9; nonlinear = true)
    @test length(unknowns(mtkcompile(sys))) == length(unknowns(mtkcompile(sys; max_linear_scc_tearing = typemax(Int))))

    # Linear in each unknown separately, but with a coefficient that contains another unknown
    sys, u = laplace_system(9; coupled = true)
    @test length(unknowns(mtkcompile(sys))) == length(unknowns(mtkcompile(sys; max_linear_scc_tearing = typemax(Int))))
end
