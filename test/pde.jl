using ModelingToolkit, DiffEqBase, LinearAlgebra, Test, Symbolics, SymbolicUtils, DomainSets

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

analytic = [u(t, x) ~ -h * x * (x - 1) * sin(x) * exp(-2 * h * t)]
analytic_function = (ps, t, x) -> -ps[1] * x * (x - 1) * sin(x) * exp(-2 * ps[1] * t)

@named pdesys = PDESystem(eq, bcs, domains, [t, x], [u(t, x)], [h => 1], analytic = analytic)
@show pdesys

@test all(isequal.(independent_variables(pdesys), [t, x]))

dx = 0:0.1:1
dt = 0:0.1:1

# Test generated analytic_func
@test all(pdesys.analytic_func[u(t, x)]([2], disct, discx) ≈
          analytic_function([2], disct, discx) for disct in dt, discx in dx)

@testset "Incompatible Array Variables" begin

    # Parameters, variables, and derivatives
    n_comp = 2
    @parameters t, x, p[1:n_comp], q[1:n_comp]
    @variables u(..)[1:n_comp]
    Dt = Differential(t)
    Dx = Differential(x)
    Dxx = Differential(x)^2
    params = reduce(vcat,[p .=> [1.5, 2.0], q .=> [1.2, 1.8]])
    # 1D PDE and boundary conditions

    eqs  = [Dt(u(t, x)[i]) ~ p[i] * Dxx(u(t, x)[i]) for i in 1:n_comp]

    bcs = [[u(0, x)[i] ~ q[i] * cos(x),
            u(t, 0)[i] ~ sin(t),
            u(t, 1)[i] ~ exp(-t) * cos(1),
            Dx(u(t,0)[i]) ~ 0.0] for i in 1:n_comp]
    bcs_collected = reduce(vcat, bcs)

    # Space and time domains
    domains = [t ∈ Interval(0.0, 1.0),
            x ∈ Interval(0.0, 1.0)]

    # PDE system
    dvs = [u(t, x)[i] for i in 1:n_comp]
    @show dvs
    @named pdesys = PDESystem(eqs, bcs_collected, domains, [t, x], dvs, params)

    # Test that the system is correctly constructed
    varname1 = Symbol("u_Any[1]")
    varname2 = Symbol("u_Any[2]")


    vars = @variables $varname1(..), $varname2(..)
    testeqs = [Dt(v(t, x)) ~ Dxx(v(t, x))*p[i] for (i, v) in enumerate(vars)]
    for i in 1:n_comp
        @test isequal(vars[i](t, x), pdesys.dvs[i])
        @test isequal(testeqs[i], pdesys.eqs[i])
    end

end
