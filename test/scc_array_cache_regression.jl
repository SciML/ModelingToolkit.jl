using ModelingToolkit, NonlinearSolve, SciMLBase, Test
using ModelingToolkitBase: t_nounits as t, D_nounits as D

scc_array_norm(x) = sum(abs2, x)
@register_symbolic scc_array_norm(x::AbstractArray)

@testset "SCC array-cache" begin
    @testset "Array across SCCs, matrix=$matrix" for matrix in (false, true)
        if matrix
            @variables m(t)[1:2, 1:1] [guess = ones(2, 1)]
        else
            @variables m(t)[1:2] [guess = ones(2)]
        end
        m1 = matrix ? m[1, 1] : m[1]
        m2 = matrix ? m[2, 1] : m[2]
        @variables z(t) [guess = 1.0]
        eqs = [
            D(m1) ~ 9 - m1^2,
            D(m2) ~ 16 - m2^2,
            z^2 ~ scc_array_norm(m),
        ]
        sys = mtkcompile(
            System(
                eqs, t; name = :array_system,
                initialization_eqs = [D(m1) ~ 0, D(m2) ~ 0]
            )
        )
        prob = InitializationProblem(sys, 0.0, Dict(); use_scc = true)
        sol = solve(prob, NewtonRaphson(); abstol = 1.0e-12, reltol = 1.0e-12)
        @test SciMLBase.successful_retcode(sol)
        @test sol[m1] ≈ 3.0
        @test sol[m2] ≈ 4.0
        @test sol[z] ≈ 5.0
    end

end
