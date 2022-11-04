using ModelingToolkit, DiffEqBase, LinearAlgebra

# Define some variables
@parameters t x
@constants h = 1
@variables u(..)
Dt = Differential(t)
Dxx = Differential(x)^2
eq = Dt(u(t, x)) ~ h * Dxx(u(t, x))
bcs = [u(0, x) ~ -h * x * (x - 1) * sin(x),
    u(t, 0) ~ 0, u(t, 1) ~ 0]

domains = [t ∈ (0.0, 1.0),
    x ∈ (0.0, 1.0)]

@named pdesys = PDESystem(eq, bcs, domains, [t, x], [u])
@show pdesys

@test all(isequal.(independent_variables(pdesys), [t, x]))
