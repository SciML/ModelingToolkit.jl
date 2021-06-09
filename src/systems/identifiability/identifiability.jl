using Symbolics, LinearAlgebra
using AbstractAlgebra
# define inputs for the example

ModelingToolkit.@parameters θ[1:4]
ModelingToolkit.@variables t, x[1:2](t), u(t), y(t)

D = ModelingToolkit.Differential(t)

n = length(x)
ℓ = length(θ)
ν = n + ℓ
ModelingToolkit.@variables Γ[1:n, 1:n](t), Λ[1:n,1:ℓ](t), W[1:n, 1:n](t)

equations = [D(x[1]) ~ x[1]^2 * θ[1] + θ[2] * x[1] * x[2] + u, D(x[2]) ~ θ[3] * x[1]^2 + θ[4] * x[1] * x[2]];
outputs = [y ~ x[1]];

## Initialize values

# random initial conditions, rational numbers
initial_conditions = Dict([x[i] => rand((1:103)) for i in 1:length(x)])

success = true
ν_cur = 1
Λ_ = zeros(Int, size(Λ))
Γ_ = Matrix{Int}(I, n, n)  

R, T = PowerSeriesRing(QQ, ν + 1, "T")
MS_n_by_n = MatrixSpace(R, n, n)
MS_n_by_ℓ = MatrixSpace(R, n, ℓ)
MS_n_by_1 = MatrixSpace(R, n, 1)


U = sum(rand(1:103) * T^j for j in 0:(ν)) # random power series

Θ = Dict([θ[i] => rand(1:102) for i in 1:length(θ)])
Φ = Dict([x[i] => R(rand(1:102)) for i in 1:length(x)])

# Building Variational System

P = [eqn.lhs - eqn.rhs  for eqn in equations];
F = [eqn.rhs  for eqn in equations]
∂P∂ẋ = ModelingToolkit.jacobian(P, D.(x))
∂P∂x = ModelingToolkit.jacobian(P, x)
∂P∂θ = ModelingToolkit.jacobian(P, θ)

∇P = [
    P,
    (∂P∂ẋ * D.(Γ)) + (∂P∂x * Γ),
    ∂P∂ẋ * D.(Λ) + ∂P∂x * (Λ) + ∂P∂θ
]

set_precision!.(values(Φ), ν_cur)
[substitute(eq, Θ) for eq in F]
A = [substitute(substitute(eq, Θ), Φ) for eq in ∂P∂ẋ]
A′ = [substitute(substitute(eq, Θ), Φ) for eq in ∂P∂x]


"""
function that solves homogeneous ODE: AW'+BW = 0
"""
function HomogeneousResolution(A, B, precision)
    # Approximate inv(A) via Newton's Method:
    
end