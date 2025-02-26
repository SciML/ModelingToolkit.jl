using ModelingToolkit, DiffEqBase, LinearAlgebra, Test
using ModelingToolkit: t_nounits as t, D_nounits as Dt

# Define some variables
@parameters x
@constants h = 1
@variables u(..)
Dxx = Differential(x)^2
eq = Dt(u(t, x)) ~ h * Dxx(u(t, x))
bcs = [u(0, x) ~ -h * x * (x - 1) * sin(x),
    u(t, 0) ~ 0, u(t, 1) ~ 0]

domains = [t ∈ (0.0, 1.0),
    x ∈ (0.0, 1.0)]

analytic = [u(t, x) ~ -h * x * (x - 1) * sin(x) * exp(-2 * h * t)]
analytic_function = (ps, t, x) -> -ps[1] * x * (x - 1) * sin(x) * exp(-2 * ps[1] * t)

@named pdesys = PDESystem(eq, bcs, domains, [t, x], [u], [h], analytic = analytic)
@show pdesys

@test all(isequal.(independent_variables(pdesys), [t, x]))

dx = 0:0.1:1
dt = 0:0.1:1

# Test generated analytic_func
@test all(pdesys.analytic_func[u(t, x)]([2], disct, discx) ≈
          analytic_function([2], disct, discx) for disct in dt, discx in dx)
