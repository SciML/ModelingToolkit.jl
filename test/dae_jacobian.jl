using ModelingToolkit
using Sundials, Test, SparseArrays
using ModelingToolkit: t_nounits as t, D_nounits as D

# Comparing solution obtained by defining explicit Jacobian function with solution obtained from
# symbolically generated Jacobian

function testjac(res, du, u, p, t) #System of equations
    res[1] = du[1] - 1.5 * u[1] + 1.0 * u[1] * u[2]
    res[2] = du[2] + 3 * u[2] - u[1] * u[2]
end

function testjac_jac(J, du, u, p, gamma, t) #Explicit Jacobian
    J[1, 1] = gamma - 1.5 + 1.0 * u[2]
    J[1, 2] = 1.0 * u[1]
    J[2, 1] = -1 * u[2]
    J[2, 2] = gamma + 3 - u[1]
    nothing
end

testjac_f = DAEFunction(testjac, jac = testjac_jac,
    jac_prototype = sparse([1, 2, 1, 2], [1, 1, 2, 2], zeros(4)))

prob1 = DAEProblem(testjac_f,
    [0.5, -2.0],
    ones(2),
    (0.0, 10.0),
    differential_vars = [true, true])
sol1 = solve(prob1, IDA(linear_solver = :KLU))

# Now MTK style solution with generated Jacobian

@variables u1(t) u2(t)
@parameters p1 p2

eqs = [D(u1) ~ p1 * u1 - u1 * u2,
    D(u2) ~ u1 * u2 - p2 * u2]

@named sys = ODESystem(eqs, t)

u0 = [u1 => 1.0,
    u2 => 1.0]

tspan = (0.0, 10.0)

du0 = [0.5, -2.0]

p = [p1 => 1.5,
    p2 => 3.0]

prob = DAEProblem(complete(sys), du0, u0, tspan, p, jac = true, sparse = true)
sol = solve(prob, IDA(linear_solver = :KLU))

@test maximum(sol - sol1) < 1e-12
