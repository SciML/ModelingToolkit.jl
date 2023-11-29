using ModelingToolkit, LinearAlgebra, SparseArrays
using Symbolics
using Symbolics: scalarize

# Define the constants for the PDE
const α₂ = 1.0
const α₃ = 1.0
const β₁ = 1.0
const β₂ = 1.0
const β₃ = 1.0
const r₁ = 1.0
const r₂ = 1.0
const _DD = 100.0
const γ₁ = 0.1
const γ₂ = 0.1
const γ₃ = 0.1
const N = 8
const X = reshape([i for i in 1:N for j in 1:N], N, N)
const Y = reshape([j for i in 1:N for j in 1:N], N, N)
const α₁ = 1.0 .* (X .>= 4 * N / 5)

const Mx = Tridiagonal([1.0 for i in 1:(N - 1)], [-2.0 for i in 1:N],
    [1.0 for i in 1:(N - 1)])
const My = copy(Mx)
Mx[2, 1] = 2.0
Mx[end - 1, end] = 2.0
My[1, 2] = 2.0
My[end, end - 1] = 2.0

# Define the initial condition as normal arrays
@variables du[1:N, 1:N, 1:3] u[1:N, 1:N, 1:3] MyA[1:N, 1:N] AMx[1:N, 1:N] DA[1:N, 1:N]

du, u, MyA, AMx, DA = scalarize.((du, u, MyA, AMx, DA))
@show typeof.((du, u, MyA, AMx, DA))

# Define the discretized PDE as an ODE function
function f(du, u, p, t)
    A = @view u[:, :, 1]
    B = @view u[:, :, 2]
    C = @view u[:, :, 3]
    dA = @view du[:, :, 1]
    dB = @view du[:, :, 2]
    dC = @view du[:, :, 3]
    mul!(MyA, My, A)
    mul!(AMx, A, Mx)
    @. DA = _DD * (MyA + AMx)
    @. dA = DA + α₁ - β₁ * A - r₁ * A * B + r₂ * C
    @. dB = α₂ - β₂ * B - r₁ * A * B + r₂ * C
    @. dC = α₃ - β₃ * C + r₁ * A * B - r₂ * C
end

f(du, u, nothing, 0.0)

multithreadedf = eval(ModelingToolkit.build_function(du, u, fillzeros = true,
    parallel = ModelingToolkit.MultithreadedForm())[2])

MyA = zeros(N, N);
AMx = zeros(N, N);
DA = zeros(N, N);
# Loop to catch synchronization issues
for i in 1:100
    _du = rand(N, N, 3)
    _u = rand(N, N, 3)
    multithreadedf(_du, _u)
    _du2 = copy(_du)
    f(_du2, _u, nothing, 0.0)
    @test _du ≈ _du2
end

#=
jac = sparse(ModelingToolkit.jacobian(vec(du),vec(u)))
fjac = eval(ModelingToolkit.build_function(jac,u,parallel=ModelingToolkit.SerialForm())[2])
multithreadedfjac = eval(ModelingToolkit.build_function(jac,u,parallel=ModelingToolkit.MultithreadedForm())[2])

u = rand(N,N,3)
J = similar(jac,Float64)
fjac(J,u)

J2 = similar(jac,Float64)
multithreadedfjac(J2,u)
@test J ≈ J2

using FiniteDiff
J3 = Array(similar(jac,Float64))
FiniteDiff.finite_difference_jacobian!(J2,(du,u)->f!(du,u,nothing,nothing),u)
maximum(J2 .- Array(J)) < 1e-5
=#

jac = ModelingToolkit.sparsejacobian(vec(du), vec(u))
serialjac = eval(ModelingToolkit.build_function(vec(jac), u)[2])
#multithreadedjac = eval(ModelingToolkit.build_function(vec(jac), u,
#    parallel = ModelingToolkit.MultithreadedForm())[2])

MyA = zeros(N, N)
AMx = zeros(N, N)
DA = zeros(N, N)
_du = rand(N, N, 3)
_u = rand(N, N, 3)

f(_du, _u, nothing, 0.0)
multithreadedf(_du, _u)

#=
using BenchmarkTools
@btime f(_du,_u,nothing,0.0)
@btime multithreadedf(_du,_u)

_jac = similar(jac,Float64)
@btime serialjac(_jac,_u)
@btime multithreadedjac(_jac,_u)
=#
