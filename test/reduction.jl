using ModelingToolkit, OrdinaryDiffEq, Test

@parameters t σ ρ β
@variables x(t) y(t) z(t) a(t) u(t) F(t)
@derivatives D'~t

test_equal(a, b) = @test isequal(simplify(a, polynorm=true), simplify(b, polynorm=true))

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y + β,
       0 ~ sin(z) - x + y,
       sin(u) ~ x + y,
       2β ~ 2,
       x ~ a,
      ]

lorenz1 = ODESystem(eqs,t,[u,x,y,z,a],[σ,ρ,β],name=:lorenz1)

lorenz1_aliased = alias_elimination(lorenz1)
reduced_eqs = [
               D(x) ~ σ * (y - x),
               D(y) ~ x*(ρ-z)-y + 1,
               0 ~ sin(z) - x + y,
               sin(u) ~ x + y,
              ]
test_equal.(equations(lorenz1_aliased), reduced_eqs)
test_equal.(states(lorenz1_aliased), [u, x, y, z])
test_equal.(observed(lorenz1_aliased), [
    β ~ 1,
    a ~ x,
])

# Multi-System Reduction

eqs1 = [
        D(x) ~ σ*(y-x) + F,
        D(y) ~ x*(ρ-z)-u,
        D(z) ~ x*y - β*z,
        u ~ x + y - z,
       ]

lorenz1 = ODESystem(eqs1,pins=[F],name=:lorenz1)

eqs2 = [
        D(x) ~ F,
        D(y) ~ x*(ρ-z)-x,
        D(z) ~ x*y - β*z,
        u ~ x - y - z
       ]

lorenz2 = ODESystem(eqs2,pins=[F],name=:lorenz2)

connected = ODESystem([lorenz2.y ~ a + lorenz1.x,
                       lorenz1.F ~ lorenz2.u,
                       lorenz2.F ~ lorenz1.u],t,[a],[],systems=[lorenz1,lorenz2])

# Reduced Flattened System

flattened_system = ModelingToolkit.flatten(connected)

aliased_flattened_system = alias_elimination(flattened_system)

@test isequal(states(aliased_flattened_system), [
        a
        lorenz1.x
        lorenz1.y
        lorenz1.z
        lorenz2.x
        lorenz2.y
        lorenz2.z
       ])

@test setdiff(parameters(aliased_flattened_system), [
        lorenz1.σ
        lorenz1.ρ
        lorenz1.β
        lorenz1.F
        lorenz2.F
        lorenz2.ρ
        lorenz2.β
       ]) |> isempty

reduced_eqs = [
               lorenz2.y ~ a + lorenz1.x, # irreducible by alias elimination
               D(lorenz1.x) ~ lorenz1.σ*(lorenz1.y-lorenz1.x) + lorenz2.x - lorenz2.y - lorenz2.z,
               D(lorenz1.y) ~ lorenz1.x*(lorenz1.ρ-lorenz1.z)-(lorenz1.x + lorenz1.y - lorenz1.z),
               D(lorenz1.z) ~ lorenz1.x*lorenz1.y - lorenz1.β*lorenz1.z,
               D(lorenz2.x) ~ lorenz1.x + lorenz1.y - lorenz1.z,
               D(lorenz2.y) ~ lorenz2.x*(lorenz2.ρ-lorenz2.z)-lorenz2.x,
               D(lorenz2.z) ~ lorenz2.x*lorenz2.y - lorenz2.β*lorenz2.z
              ]
test_equal.(equations(aliased_flattened_system), reduced_eqs)

observed_eqs = [
                lorenz1.F ~ lorenz2.u,
                lorenz2.F ~ lorenz1.u,
                lorenz1.u ~ lorenz1.x + lorenz1.y - lorenz1.z,
                lorenz2.u ~ lorenz2.x - lorenz2.y - lorenz2.z,
               ]
test_equal.(observed(aliased_flattened_system), observed_eqs)

# issue #578

let
    @variables t x(t) y(t) z(t);
    @derivatives D'~t;
    eqs = [
        D(x) ~ x + y
        x ~ y
    ];
    sys = ODESystem(eqs, t, [x], []);
    asys = alias_elimination(flatten(sys))

    test_equal.(asys.eqs, [D(x) ~ 2x])
    test_equal.(asys.observed, [y ~ x])
end

# issue #716
let
    @parameters t
    @derivatives D'~t
    @variables x(t), u(t), y(t)
    @parameters a, b, c, d
    ol = ODESystem([D(x) ~ a * x + b * u, y ~ c * x], t, name=:ol)
    @variables u_c(t), y_c(t)
    @parameters k_P
    pc = ODESystem(Equation[], t, pins=[y_c], observed = [u_c ~ k_P * y_c], name=:pc)
    connections = [
                   ol.u ~ pc.u_c,
                   y_c ~ ol.y
                  ]
    connected = ODESystem(connections, t, systems=[ol, pc])

    @test equations(connected) isa Vector{Equation}
    @test_nowarn flatten(connected)
end
