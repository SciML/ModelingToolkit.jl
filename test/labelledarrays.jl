using ModelingToolkit, StaticArrays, LinearAlgebra, LabelledArrays
using DiffEqBase, ForwardDiff
using Test
using ModelingToolkit: t_nounits as t, D_nounits as D

# Define some variables
@parameters σ ρ β
@variables x(t) y(t) z(t)

# Define a differential equation
eqs = [D(x) ~ σ * (y - x),
    D(y) ~ t * x * (ρ - z) - y,
    D(z) ~ x * y - β * z]

@named de = ODESystem(eqs, t)
de = complete(de)
ff = ODEFunction(de, [x, y, z], [σ, ρ, β], jac = true)

a = @SVector [1.0, 2.0, 3.0]
b = SLVector(x = 1.0, y = 2.0, z = 3.0)
c = [1.0, 2.0, 3.0]
p = (SLVector(σ = 10.0, ρ = 26.0, β = 8 / 3),)
@test ff(a, p, 0.0) isa SVector
@test typeof(ff(b, p, 0.0)) <: SLArray
@test ff(c, p, 0.0) isa Vector
@test ff(a, p, 0.0) == ff(b, p, 0.0)
@test ff(a, p, 0.0) == ff(c, p, 0.0)

@test ff.jac(a, p, 0.0) isa SMatrix
@test typeof(ff.jac(b, p, 0.0)) <: SMatrix
@test ff.jac(c, p, 0.0) isa Matrix
@test ff.jac(a, p, 0.0) == ff.jac(b, p, 0.0)
@test ff.jac(a, p, 0.0) == ff.jac(c, p, 0.0)

# Test similar_type
@test ff(b, p, ForwardDiff.Dual(0.0, 1.0)) isa SLArray
d = LVector(x = 1.0, y = 2.0, z = 3.0)
@test ff(d, p, ForwardDiff.Dual(0.0, 1.0)) isa LArray
@test ff.jac(b, p, ForwardDiff.Dual(0.0, 1.0)) isa SArray
@test eltype(ff.jac(b, p, ForwardDiff.Dual(0.0, 1.0))) <: ForwardDiff.Dual
@test ff.jac(d, p, ForwardDiff.Dual(0.0, 1.0)) isa Array
@inferred ff.jac(d, p, ForwardDiff.Dual(0.0, 1.0))
@test eltype(ff.jac(d, p, ForwardDiff.Dual(0.0, 1.0))) <: ForwardDiff.Dual

## https://github.com/SciML/ModelingToolkit.jl/issues/1054
using LabelledArrays
using ModelingToolkit

# ODE model: simple SIR model with seasonally forced contact rate
function SIR!(du, u, p, t)

    # Unknowns
    (S, I, R) = u[1:3]
    N = S + I + R

    # params
    β = p.β
    η = p.η
    φ = p.φ
    ω = 1.0 / p.ω
    μ = p.μ
    σ = p.σ

    # FOI
    βeff = β * (1.0 + η * cos(2.0 * π * (t - φ) / 365.0))
    λ = βeff * I / N

    # change in unknowns
    du[1] = (μ * N - λ * S - μ * S + ω * R)
    du[2] = (λ * S - σ * I - μ * I)
    du[3] = (σ * I - μ * R - ω * R)
    du[4] = (σ * I) # cumulative incidence
end

# Solver settings
tmin = 0.0
tmax = 10.0 * 365.0
tspan = (tmin, tmax)

# Initiate ODE problem
theta_fix = [1.0 / (80 * 365)]
theta_est = [0.28, 0.07, 1.0 / 365.0, 1.0, 1.0 / 5.0]
p = @LArray [theta_est; theta_fix] (:β, :η, :ω, :φ, :σ, :μ)
u0 = @LArray [9998.0, 1.0, 1.0, 1.0] (:S, :I, :R, :C)

# Initiate ODE problem
problem = ODEProblem(SIR!, u0, tspan, p)
sys = complete(modelingtoolkitize(problem))

@test all(any(isequal(x), parameters(sys))
for x in ModelingToolkit.unwrap.(@variables(β, η, ω, φ, σ, μ)))
@test all(isequal.(Symbol.(unknowns(sys)), Symbol.(@variables(S(t), I(t), R(t), C(t)))))
