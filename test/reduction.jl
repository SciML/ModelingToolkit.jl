using ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkit: topsort_equations

@variables t x(t) y(t) z(t) k(t)
eqs = [
       x ~ y + z
       z ~ 2
       y ~ 2z + k
      ]

sorted_eq = topsort_equations(eqs, [x, y, z, k])

ref_eq = [
          z ~ 2
          y ~ 2z + k
          x ~ y + z
         ]
@test ref_eq == sorted_eq

@test_throws ArgumentError topsort_equations([
                                             x ~ y + z
                                             z ~ 2
                                             y ~ 2z + x
                                            ], [x, y, z, k])

@parameters t σ ρ β
@variables x(t) y(t) z(t) a(t) u(t) F(t)
D = Differential(t)

test_equal(a, b) = @test isequal(simplify(a, polynorm=true), simplify(b, polynorm=true))

eqs = [
       D(x) ~ σ*(y-x)
       D(y) ~ x*(ρ-z)-y + β
       0 ~ sin(z) - x + y
       sin(u) ~ x + y
       2β ~ 2
       x ~ a
      ]

lorenz1 = ODESystem(eqs,t,name=:lorenz1)

lorenz1_aliased = alias_elimination(lorenz1)
reduced_eqs = [
               D(x) ~ σ * (y - x),
               D(y) ~ x*(ρ-z)-y + β,
               0 ~ sin(z) - x + y,
               0 ~ x + y - sin(u),
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

lorenz = name -> ODESystem(eqs1,t,pins=[F],name=name)
lorenz1 = lorenz(:lorenz1)
lorenz2 = lorenz(:lorenz2)

connected = ODESystem([s ~ a + lorenz1.x
                       lorenz2.y ~ s
                       lorenz1.F ~ lorenz2.u
                       lorenz2.F ~ lorenz1.u],t,systems=[lorenz1,lorenz2])

# Reduced Flattened System

flattened_system = ModelingToolkit.flatten(connected)

aliased_flattened_system = alias_elimination(flattened_system)

@test setdiff(states(aliased_flattened_system), [
        a
        lorenz1.x
        lorenz1.y
        lorenz1.z
        lorenz2.x
        lorenz2.y
        lorenz2.z
       ]) |> isempty

@test setdiff(parameters(aliased_flattened_system), [
        lorenz1.σ
        lorenz1.ρ
        lorenz1.β
        lorenz2.σ
        lorenz2.ρ
        lorenz2.β
       ]) |> isempty

reduced_eqs = [
               0 ~ s - lorenz2.y
               D(lorenz1.x) ~ lorenz1.F + lorenz1.σ*(lorenz1.y + -1lorenz1.x)
               D(lorenz1.y) ~ -1lorenz1.u + lorenz1.x*(lorenz1.ρ + -1lorenz1.z)
               D(lorenz1.z) ~ lorenz1.x*lorenz1.y + -1lorenz1.β*lorenz1.z
               D(lorenz2.x) ~ lorenz2.F + lorenz2.σ*(lorenz2.y + -1lorenz2.x)
               D(lorenz2.y) ~ -1lorenz2.u + lorenz2.x*(lorenz2.ρ + -1lorenz2.z)
               D(lorenz2.z) ~ lorenz2.x*lorenz2.y + -1lorenz2.β*lorenz2.z
              ]
test_equal.(equations(aliased_flattened_system), reduced_eqs)

observed_eqs = [
                s ~ a + lorenz1.x
                lorenz1.u ~ lorenz1.x + lorenz1.y - lorenz1.z
                lorenz2.u ~ lorenz2.x + lorenz2.y - lorenz2.z
                lorenz2.F ~ lorenz1.u
                lorenz1.F ~ lorenz2.u
               ]
test_equal.(observed(aliased_flattened_system), observed_eqs)

# issue #578

let
    @variables t x(t) y(t) z(t)
    D = Differential(t)
    eqs = [
        D(x) ~ x + y
        x ~ y
    ]
    sys = ODESystem(eqs, t)
    asys = alias_elimination(flatten(sys))

    test_equal.(asys.eqs, [D(x) ~ x + y])
    test_equal.(asys.observed, [y ~ x])
end

# issue #724 and #716
let
    @parameters t
    D = Differential(t)
    @variables x(t), u(t), y(t)
    @parameters a, b, c, d
    ol = ODESystem([D(x) ~ a * x + b * u; y ~ c * x + d * u], t, pins=[u], name=:ol)
    @variables u_c(t), y_c(t)
    @parameters k_P
    pc = ODESystem(Equation[u_c ~ k_P * y_c], t, pins=[y_c], name=:pc)
    connections = [
        ol.u ~ pc.u_c
        pc.y_c ~ ol.y
    ]
    connected = ODESystem(connections, t, systems=[ol, pc])
    @test equations(connected) isa Vector{Equation}
    sys = flatten(connected)
    reduced_sys = alias_elimination(sys)
    ref_eqs = [
               D(ol.x) ~ ol.a*ol.x + ol.b*ol.u
               0 ~ ol.c*ol.x + ol.d*ol.u + -1ol.y
               0 ~ pc.k_P*pc.y_c + -1pc.u_c
              ]
    @test ref_eqs == equations(reduced_sys)
end
