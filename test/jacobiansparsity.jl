using ModelingToolkit, SparseArrays, OrdinaryDiffEq, DiffEqBase

N = 3
xyd_brusselator = range(0, stop = 1, length = N)
brusselator_f(x, y, t) = (((x - 0.3)^2 + (y - 0.6)^2) <= 0.1^2) * (t >= 1.1) * 5.0
lim(a, N) = ModelingToolkit.ifelse(a == N + 1, 1, ModelingToolkit.ifelse(a == 0, N, a))
function brusselator_2d_loop(du, u, p, t)
    A, B, alpha, dx = p
    alpha = alpha / dx^2
    @inbounds for I in CartesianIndices((N, N))
        i, j = Tuple(I)
        x, y = xyd_brusselator[I[1]], xyd_brusselator[I[2]]
        ip1, im1, jp1, jm1 = lim(i + 1, N), lim(i - 1, N), lim(j + 1, N),
        lim(j - 1, N)
        du[i,
        j,
        1] = alpha * (u[im1, j, 1] + u[ip1, j, 1] + u[i, jp1, 1] + u[i, jm1, 1] -
                       4u[i, j, 1]) +
                      B + u[i, j, 1]^2 * u[i, j, 2] - (A + 1) * u[i, j, 1] +
                      brusselator_f(x, y, t)
        du[i,
        j,
        2] = alpha * (u[im1, j, 2] + u[ip1, j, 2] + u[i, jp1, 2] + u[i, jm1, 2] -
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
#@test_nowarn solve(prob, Rosenbrock23())
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

@testset "W matrix sparsity" begin
    t = ModelingToolkit.t_nounits
    D = ModelingToolkit.D_nounits
    @parameters g
    @variables x(t) y(t) λ(t)
    eqs = [D(D(x)) ~ λ * x
           D(D(y)) ~ λ * y - g
           x^2 + y^2 ~ 1]
    @mtkcompile pend = System(eqs, t)

    u0 = [x => 1, y => 0]
    prob = ODEProblem(
        pend, [u0; [g => 1]], (0, 11.5), guesses = [λ => 1], sparse = true, jac = true)
    jac, jac! = generate_jacobian(pend; expression = Val{false}, sparse = true)
    jac_prototype = ModelingToolkit.jacobian_sparsity(pend)
    W_prototype = ModelingToolkit.W_sparsity(pend)
    @test nnz(W_prototype) == nnz(jac_prototype) + 2

    # jac_prototype should be the same as W_prototype
    @test findnz(prob.f.jac_prototype)[1:2] == findnz(W_prototype)[1:2]

    u = zeros(5)
    p = prob.p
    t = 0.0
    @test_throws AssertionError jac!(similar(jac_prototype, Float64), u, p, t)

    W, W! = generate_W(pend; expression = Val{false}, sparse = true)
    γ = 0.1
    M = sparse(calculate_massmatrix(pend))
    @test_throws AssertionError W!(similar(jac_prototype, Float64), u, p, γ, t)
    @test W!(similar(W_prototype, Float64), u, p, γ, t) ==
          0.1 * M + jac!(similar(W_prototype, Float64), u, p, t)
end

@testset "Issue#3556: Numerical accuracy" begin
    t = ModelingToolkit.t_nounits
    D = ModelingToolkit.D_nounits
    @parameters g
    @variables x(t) y(t) [state_priority = 10] λ(t)
    eqs = [D(D(x)) ~ λ * x
           D(D(y)) ~ λ * y - g
           x^2 + y^2 ~ 1]
    @mtkcompile pend = System(eqs, t)
    prob = ODEProblem(pend, [x => 0.0, D(x) => 1.0, g => 1.0], (0.0, 1.0);
        guesses = [y => 1.0, λ => 1.0], jac = true, sparse = true)
    J = deepcopy(prob.f.jac_prototype)
    prob.f.jac(J, prob.u0, prob.p, 1.0)
    # this currently works but may not continue to do so
    # see https://github.com/SciML/ModelingToolkit.jl/pull/3556#issuecomment-2792664039
    @test J == prob.f.jac(prob.u0, prob.p, 1.0)
    @test J ≈ prob.f.jac(prob.u0, prob.p, 1.0)
end

# https://github.com/SciML/ModelingToolkit.jl/issues/3871
@testset "Issue#3871: Sparsity with observed derivatives" begin
    t = ModelingToolkit.t_nounits
    D = ModelingToolkit.D_nounits
    @variables x(t) y(t)
    @mtkcompile sys = System([D(x) ~ x * D(y), D(y) ~ x - y], t)
    @test ModelingToolkit.jacobian_sparsity(sys) == [1 1; 1 1] # all nonzero
    J1 = calculate_jacobian(sys)
    J2 = isequal(unknowns(sys)[1], x) ? [2x-y -x; 1 -1] : [-1 1; -x 2x-y] # analytical result
    @test isequal(J1, J2)
    prob = ODEProblem(sys, [x => 1.0, y => 0.0], (0.0, 1.0); jac = true, sparse = true)
    sol = solve(prob, FBDF())
    @test SciMLBase.successful_retcode(sol)
    ts = ModelingToolkit.get_tearing_state(sys)
    for ieq in 1:2
        vars1 = ts.fullvars[ModelingToolkit.BipartiteGraphs.𝑠neighbors(ts.structure.graph, ieq)]
        vars2 = ModelingToolkit.vars(equations(sys)[ieq])
        @test issetequal(vars1, vars2)
    end
end
