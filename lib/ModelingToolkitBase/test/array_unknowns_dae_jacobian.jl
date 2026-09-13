using ModelingToolkitBase, SciMLBase, Test, SparseArrays

@parameters t
D = Differential(t)

@testset "DAE Jacobians with array unknowns" begin
    for shape in ((3,), (2, 2)), array_unknowns in (false, true),
            array_equations in (false, true)
        @testset "shape=$shape array_unknowns=$array_unknowns array_equations=$array_equations" begin
            @variables x(t)[(1:s for s in shape)...] z(t)
            xs = vec(collect(x))
            if array_equations
                eqs = [D(x) ~ -2x .+ z, z ~ sum(x .^ 2)]
                vars = array_unknowns ? [x, z] : vcat(xs, z)
                n = length(xs) + 1
            else
                n = length(xs)
                eqs = [D(xs[i]) ~ -2xs[i] + xs[i + 1] for i in 1:(n - 1)]
                push!(eqs, xs[n] ~ xs[1]^2)
                vars = array_unknowns ? [x] : xs
            end
            @named raw = System(eqs, t, vars, [])
            sys = complete(raw)
            u = Float64.(1:n)
            du = fill(0.5, n)
            gamma = 0.7
            expected = zeros(n, n)
            for i in 1:(n - 1)
                expected[i, i] = -2 - gamma
                expected[i, array_equations ? n : i + 1] = 1
            end
            expected[n, n] = -1
            if array_equations
                expected[n, 1:(n - 1)] .= 2u[1:(n - 1)]
                expected_residual = vcat(
                    -2u[1:(n - 1)] .+ u[n] .- du[1:(n - 1)],
                    sum(abs2, u[1:(n - 1)]) - u[n]
                )
            else
                expected[n, 1] = 2u[1]
                expected_residual = vcat(
                    -2u[1:(n - 1)] + u[2:n] - du[1:(n - 1)], u[1]^2 - u[n]
                )
            end
            for analytic in (false, true), sparse in (false, true)
                @testset "jac=$analytic sparse=$sparse" begin
                    f = DAEFunction(sys; jac = analytic, sparse, u0 = u)
                    residual = zeros(n)
                    f(residual, du, u, nothing, 0.0)
                    @test residual ≈ expected_residual
                    if sparse
                        @test size(f.jac_prototype) == (n, n)
                        @test findnz(f.jac_prototype)[1:2] ==
                            findnz(SparseArrays.sparse(expected))[1:2]
                    end
                    if analytic
                        out = sparse ? copy(f.jac_prototype) : zeros(n, n)
                        f.jac(out, du, u, nothing, gamma, 0.0)
                        @test Matrix(out) ≈ expected
                        @test Matrix(f.jac(du, u, nothing, gamma, 0.0)) ≈ expected
                    end
                end
            end
        end
    end
end
