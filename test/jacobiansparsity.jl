using OrdinaryDiffEq, ModelingToolkit, Test, SparseArrays

N = 3
xyd_brusselator = range(0, stop = 1, length = N)
brusselator_f(x, y, t) = (((x - 0.3)^2 + (y - 0.6)^2) <= 0.1^2) * (t >= 1.1) * 5.0
limit(a, N) = ModelingToolkit.ifelse(a == N + 1, 1, ModelingToolkit.ifelse(a == 0, N, a))
function brusselator_2d_loop(du, u, p, t)
    A, B, alpha, dx = p
    alpha = alpha / dx^2
    @inbounds for I in CartesianIndices((N, N))
        i, j = Tuple(I)
        x, y = xyd_brusselator[I[1]], xyd_brusselator[I[2]]
        ip1, im1, jp1, jm1 = limit(i + 1, N), limit(i - 1, N), limit(j + 1, N),
        limit(j - 1, N)
        du[i, j, 1] = alpha * (u[im1, j, 1] + u[ip1, j, 1] + u[i, jp1, 1] + u[i, jm1, 1] -
                       4u[i, j, 1]) +
                      B + u[i, j, 1]^2 * u[i, j, 2] - (A + 1) * u[i, j, 1] +
                      brusselator_f(x, y, t)
        du[i, j, 2] = alpha * (u[im1, j, 2] + u[ip1, j, 2] + u[i, jp1, 2] + u[i, jm1, 2] -
                       4u[i, j, 2]) +
                      A * u[i, j, 1] - u[i, j, 1]^2 * u[i, j, 2]
    end
end

# Test with tuple parameters
p = (3.4, 1.0, 10.0, step(xyd_brusselator))

function init_brusselator_2d(xyd)
    N = length(xyd)
    u = zeros(N, N, 2)
    for I in CartesianIndices((N, N))
        x = xyd[I[1]]
        y = xyd[I[2]]
        u[I, 1] = 22 * (y * (1 - y))^(3 / 2)
        u[I, 2] = 27 * (x * (1 - x))^(3 / 2)
    end
    u
end

u0 = init_brusselator_2d(xyd_brusselator)
prob_ode_brusselator_2d = ODEProblem(brusselator_2d_loop,
    u0, (0.0, 11.5), p)
sys = complete(modelingtoolkitize(prob_ode_brusselator_2d))

# test sparse jacobian pattern only.
prob = ODEProblem(sys, u0, (0, 11.5), sparse = true, jac = false)
JP = prob.f.jac_prototype
@test findnz(Symbolics.jacobian_sparsity(map(x -> x.rhs, equations(sys)),
    unknowns(sys)))[1:2] ==
      findnz(JP)[1:2]

# test sparse jacobian
prob = ODEProblem(sys, u0, (0, 11.5), sparse = true, jac = true)
@test_nowarn solve(prob, Rosenbrock23())
@test findnz(calculate_jacobian(sys, sparse = true))[1:2] ==
      findnz(prob.f.jac_prototype)[1:2]

# test when not sparse
prob = ODEProblem(sys, u0, (0, 11.5), sparse = false, jac = true)
@test prob.f.jac_prototype == nothing

prob = ODEProblem(sys, u0, (0, 11.5), sparse = false, jac = false)
@test prob.f.jac_prototype == nothing

# test when u0 is nothing
f = DiffEqBase.ODEFunction(sys, u0 = nothing, sparse = true, jac = true)
@test findnz(f.jac_prototype)[1:2] == findnz(JP)[1:2]
@test eltype(f.jac_prototype) == Float64

f = DiffEqBase.ODEFunction(sys, u0 = nothing, sparse = true, jac = false)
@test findnz(f.jac_prototype)[1:2] == findnz(JP)[1:2]
@test eltype(f.jac_prototype) == Float64

# test when u0 is not Float64
u0 = similar(init_brusselator_2d(xyd_brusselator), Float32)
prob_ode_brusselator_2d = ODEProblem(brusselator_2d_loop,
    u0, (0.0, 11.5), p)
sys = complete(modelingtoolkitize(prob_ode_brusselator_2d))

prob = ODEProblem(sys, u0, (0, 11.5), sparse = true, jac = false)
@test eltype(prob.f.jac_prototype) == Float32

prob = ODEProblem(sys, u0, (0, 11.5), sparse = true, jac = true)
@test eltype(prob.f.jac_prototype) == Float32
