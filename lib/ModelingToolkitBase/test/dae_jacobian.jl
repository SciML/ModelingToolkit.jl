using ModelingToolkitBase
using Sundials, Test, SparseArrays
using ModelingToolkitBase: t_nounits as t, D_nounits as D

# Comparing solution obtained by defining explicit Jacobian function with solution obtained from
# symbolically generated Jacobian

function testjac(res, du, u, p, t) #System of equations
    res[1] = du[1] - 1.5 * u[1] + 1.0 * u[1] * u[2]
    return res[2] = du[2] + 3 * u[2] - u[1] * u[2]
end

function testjac_jac(J, du, u, p, gamma, t) #Explicit Jacobian
    J[1, 1] = gamma - 1.5 + 1.0 * u[2]
    J[1, 2] = 1.0 * u[1]
    J[2, 1] = -1 * u[2]
    J[2, 2] = gamma + 3 - u[1]
    return nothing
end

testjac_f = DAEFunction(
    testjac, jac = testjac_jac,
    jac_prototype = sparse([1, 2, 1, 2], [1, 1, 2, 2], zeros(4))
)

prob1 = DAEProblem(
    testjac_f,
    [0.5, -2.0],
    ones(2),
    (0.0, 10.0),
    differential_vars = [true, true]
)
sol1 = solve(prob1, IDA(linear_solver = :KLU))

# Now MTK style solution with generated Jacobian

@variables u1(t) u2(t)
@parameters p1 p2

eqs = [
    D(u1) ~ p1 * u1 - u1 * u2,
    D(u2) ~ u1 * u2 - p2 * u2,
]

@named sys = System(eqs, t)

u0 = [
    u1 => 1.0,
    u2 => 1.0,
]

tspan = (0.0, 10.0)

du0 = [D(u1) => 0.5, D(u2) => -2.0]

p = [
    p1 => 1.5,
    p2 => 3.0,
]

prob = DAEProblem(complete(sys), [du0; p], tspan, jac = true, sparse = true, missing_guess_value = MissingGuessValue.Constant(1.0))
sol = solve(prob, IDA(linear_solver = :KLU))

@test maximum(sol - sol1) < 2.0e-12

@testset "DAE Jacobians include both sides of scalar equations" begin
    @variables x(t) y(t)
    for differential_eq in (D(x) ~ y, 0 ~ y - D(x)),
            algebraic_eq in (y ~ x^2, 0 ~ x^2 - y)
        @named raw = System([differential_eq, algebraic_eq], t, [x, y], [])
        sys = complete(raw)
        u, du = [2.0, 5.0], [1.0, 0.0]
        for analytic in (false, true), sparse in (false, true)
            @testset "jac=$analytic sparse=$sparse" begin
                f = DAEFunction(sys; jac = analytic, sparse, u0 = u)
                residual = zeros(2)
                f(residual, du, u, nothing, 0.0)
                @test residual == [4.0, -1.0]
                if sparse
                    @test size(f.jac_prototype) == (2, 2)
                    @test findnz(f.jac_prototype)[1:2] == ([1, 2, 1, 2], [1, 1, 2, 2])
                end
                if analytic
                    for gamma in (0.0, 0.7, 1.3)
                        expected = [-gamma 1.0; 2u[1] -1.0]
                        jac = sparse ? copy(f.jac_prototype) : zeros(2, 2)
                        f.jac(jac, du, u, nothing, gamma, 0.0)
                        @test Matrix(jac) ≈ expected
                        @test Matrix(f.jac(du, u, nothing, gamma, 0.0)) ≈ expected
                    end
                end
            end
        end
    end
end
