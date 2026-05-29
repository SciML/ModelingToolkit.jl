using ModelingToolkitBase, DiffEqBase, LinearAlgebra, Test
using ModelingToolkitBase: t_nounits as t, D_nounits as Dt
using Symbolics: getname

# Define some variables
@parameters x
@constants h = 1
@variables u(..)
Dxx = Differential(x)^2
eq = Dt(u(t, x)) ~ h * Dxx(u(t, x))
bcs = [
    u(0, x) ~ -h * x * (x - 1) * sin(x),
    u(t, 0) ~ 0, u(t, 1) ~ 0,
]

domains = [
    t ∈ (0.0, 1.0),
    x ∈ (0.0, 1.0),
]

analytic = [u(t, x) ~ -h * x * (x - 1) * sin(x) * exp(-2 * h * t)]
analytic_function = (ps, t, x) -> -ps[1] * x * (x - 1) * sin(x) * exp(-2 * ps[1] * t)

@named pdesys = PDESystem(eq, bcs, domains, [t, x], [u], [h], analytic = analytic)
@show pdesys

@test all(isequal.(independent_variables(pdesys), [t, x]))

dx = 0:0.1:1
dt = 0:0.1:1

# Test generated analytic_func
@test all(
    pdesys.analytic_func[u(t, x)]([2], disct, discx) ≈
        analytic_function([2], disct, discx) for disct in dt, discx in dx
)

# Test sys.x accessor pattern for PDESystem
using ModelingToolkitBase: renamespace, does_namespacing

@testset "PDESystem property accessor (sys.x)" begin
    @parameters y α β
    @variables v(..) w(..)
    Dy = Differential(y)
    eq2 = [
        Dt(v(t, y)) ~ α * Dy(Dy(v(t, y))),
        Dt(w(t, y)) ~ β * Dy(Dy(w(t, y))),
    ]
    bcs2 = [
        v(0, y) ~ sin(y), w(0, y) ~ cos(y),
        v(t, 0) ~ 0, v(t, 1) ~ 0,
        w(t, 0) ~ 1, w(t, 1) ~ 1,
    ]
    domains2 = [t ∈ (0.0, 1.0), y ∈ (0.0, 1.0)]

    @named pdesys2 = PDESystem(eq2, bcs2, domains2, [t, y], [v, w], [α, β])

    # PDESystem does namespacing by default (not complete), matching System behavior
    @test does_namespacing(pdesys2)

    # Access dependent variables (namespaced, matching System behavior)
    @test !isequal(v, pdesys2.v)
    @test isequal(renamespace(pdesys2, v), pdesys2.v)
    @test !isequal(w, pdesys2.w)
    @test isequal(renamespace(pdesys2, w), pdesys2.w)

    # Access parameters (namespaced)
    @test !isequal(α, pdesys2.α)
    @test isequal(renamespace(pdesys2, α), pdesys2.α)
    @test !isequal(β, pdesys2.β)
    @test isequal(renamespace(pdesys2, β), pdesys2.β)

    # Access independent variables — IVs are not namespaced
    @test isequal(t, pdesys2.t)
    @test isequal(y, pdesys2.y)

    # Nonexistent variable throws
    @test_throws ArgumentError pdesys2.nonexistent

    # propertynames lists all accessible names
    pnames = propertynames(pdesys2)
    @test :v ∈ pnames
    @test :w ∈ pnames
    @test :α ∈ pnames
    @test :β ∈ pnames
    @test :t ∈ pnames
    @test :y ∈ pnames

    # Struct field access still works
    @test pdesys2.eqs isa Vector
    @test length(pdesys2.eqs) == 2
    @test pdesys2.bcs isa Vector
end

@testset "PDESystem property accessor without parameters" begin
    @parameters z
    @variables f(..)
    Dz = Differential(z)
    eq3 = Dt(f(t, z)) ~ Dz(Dz(f(t, z)))
    bcs3 = [f(0, z) ~ sin(z), f(t, 0) ~ 0, f(t, 1) ~ 0]
    domains3 = [t ∈ (0.0, 1.0), z ∈ (0.0, 1.0)]

    @named pdesys3 = PDESystem(eq3, bcs3, domains3, [t, z], [f])

    # Access dependent variable (namespaced)
    @test isequal(renamespace(pdesys3, f), pdesys3.f)

    # Access independent variables (not namespaced)
    @test isequal(t, pdesys3.t)
    @test isequal(z, pdesys3.z)
end

@testset "PDESystem property accessor with subsystems" begin
    @parameters s2 γ
    @variables p(..) q(..)
    Ds2 = Differential(s2)

    eq_inner = [Dt(p(t, s2)) ~ γ * Ds2(p(t, s2))]
    bcs_inner = [p(0, s2) ~ sin(s2), p(t, 0) ~ 0, p(t, 1) ~ 0]
    domains_inner = [t ∈ (0.0, 1.0), s2 ∈ (0.0, 1.0)]
    @named inner = PDESystem(eq_inner, bcs_inner, domains_inner, [t, s2], [p], [γ])

    eq_outer = [Dt(q(t, s2)) ~ Ds2(q(t, s2))]
    bcs_outer = [q(0, s2) ~ cos(s2), q(t, 0) ~ 1, q(t, 1) ~ 1]
    domains_outer = [t ∈ (0.0, 1.0), s2 ∈ (0.0, 1.0)]
    @named outer = PDESystem(
        eq_outer, bcs_outer, domains_outer, [t, s2], [q], [γ];
        systems = [inner]
    )

    # Subsystem should be listed in propertynames
    @test :inner ∈ propertynames(outer)

    # Access own variables (namespaced)
    @test isequal(renamespace(outer, q), outer.q)
end
