using Pkg
Pkg.activate("../../../")
using ModelingToolkit
using Symbolics

@parameters t, theta[1:4]
@variables x[1:2](t), y
D = Differential(t)

n = length(x)
ℓ = length(θ)

ode = [D(x[1]) ~ x[1]^2*θ[1] + θ[2]*x[1]*x[2]+U, D(x[2])~θ[3]*x[1]^2+θ[4]*x[1]*x[2]]
outputs = [y~x[1]]

# How to generate power series solution?



# constructing P 
@variables ẋ[1:length(x)](t)

subs=[substitute(expr.lhs, Dict([D(x[i]) => ẋ[i] for i in 1:length(x)]))~expr.rhs for expr in ode]

P=[x.lhs-x.rhs for x in subs]

# constructing ∇P
dPdẋ = Symbolics.jacobian(P, ẋ)
dPdx = Symbolics.jacobian(P, x)
dPdθ = Symbolics.jacobian(P, θ)



