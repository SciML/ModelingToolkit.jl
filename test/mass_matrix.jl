using OrdinaryDiffEq, ModelingToolkit, Test, LinearAlgebra, StaticArrays
using ModelingToolkit: t_nounits as t, D_nounits as D, MTKParameters

@variables y(t)[1:3]
@parameters k[1:3]

eqs = [D(y[1]) ~ -k[1] * y[1] + k[3] * y[2] * y[3],
    D(y[2]) ~ k[1] * y[1] - k[3] * y[2] * y[3] - k[2] * y[2]^2,
    0 ~ y[1] + y[2] + y[3] - 1]

@named sys = System(eqs, t, collect(y), [k])
sys = complete(sys)
@test_throws ModelingToolkit.OperatorIndepvarMismatchError System(eqs, y[1])
M = calculate_massmatrix(sys)
@test M isa Diagonal
@test M == [1 0 0
            0 1 0
            0 0 0]

prob_mm = ODEProblem(sys, [y => [1.0, 0.0, 0.0], k => [0.04, 3e7, 1e4]], (0.0, 1e5))
@test prob_mm.f.mass_matrix isa Diagonal{Float64, Vector{Float64}}
sol = solve(prob_mm, Rodas5(), reltol = 1e-8, abstol = 1e-8)
prob_mm = ODEProblem(sys, SA[y => [1.0, 0.0, 0.0], k => [0.04, 3e7, 1e4]], (0.0, 1e5))
@test prob_mm.f.mass_matrix isa Diagonal{Float64, SVector{3, Float64}}

function rober(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du[3] = y₁ + y₂ + y₃ - 1
    nothing
end
f = ODEFunction(rober, mass_matrix = M)
prob_mm2 = ODEProblem(f, [1.0, 0.0, 0.0], (0.0, 1e5), (0.04, 3e7, 1e4))
sol2 = solve(prob_mm2, Rodas5(), reltol = 1e-8, abstol = 1e-8, tstops = sol.t,
    adaptive = false)

# MTK expression are canonicalized, so the floating point numbers are slightly
# different
@test Array(sol) ≈ Array(sol2)

# Test mass matrix in the identity case
eqs = [D(y[1]) ~ y[1], D(y[2]) ~ y[2], D(y[3]) ~ y[3]]

@named sys = System(eqs, t, collect(y), [k])

@test calculate_massmatrix(sys) === I

@testset "Mass matrix `isa Diagonal` for `SDEProblem`" begin
    eqs = [D(y[1]) ~ -k[1] * y[1] + k[3] * y[2] * y[3],
        D(y[2]) ~ k[1] * y[1] - k[3] * y[2] * y[3] - k[2] * y[2]^2,
        0 ~ y[1] + y[2] + y[3] - 1]

    @named sys = System(eqs, t, collect(y), [k])
    @named sys = SDESystem(sys, [1, 1, 0])
    sys = complete(sys)
    prob = SDEProblem(sys, [y => [1.0, 0.0, 0.0], k => [0.04, 3e7, 1e4]], (0.0, 1e5))
    @test prob.f.mass_matrix isa Diagonal{Float64, Vector{Float64}}
    prob = SDEProblem(sys, SA[y => [1.0, 0.0, 0.0], k => [0.04, 3e7, 1e4]], (0.0, 1e5))
    @test prob.f.mass_matrix isa Diagonal{Float64, SVector{3, Float64}}
end
