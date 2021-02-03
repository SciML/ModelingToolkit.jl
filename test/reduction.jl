using ModelingToolkit, OrdinaryDiffEq, Test, NonlinearSolve
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
       2x ~ 3a
       2u ~ 3x
      ]

lorenz1 = ODESystem(eqs,t,name=:lorenz1)

lorenz1_aliased = alias_elimination(lorenz1)
reduced_eqs = [
               D(x) ~ σ*(y - x)
               D(y) ~ β + x*(ρ - z) - y
               0 ~ y + sin(z) - x
               0 ~ x + y - sin(1.5x)
              ]
test_equal.(equations(lorenz1_aliased), reduced_eqs)
@test isempty(setdiff(states(lorenz1_aliased), [u, x, y, z]))
test_equal.(observed(lorenz1_aliased), [
                                        a ~ 2/3 * x,
                                        u ~ 3/2 * x,
                                       ])

# Multi-System Reduction

@variables s(t)
eqs1 = [
        D(x) ~ σ*(y-x) + F,
        D(y) ~ x*(ρ-z)-u,
        D(z) ~ x*y - β*z,
        u ~ x + y - z,
       ]

lorenz = name -> ODESystem(eqs1,t,name=name)
lorenz1 = lorenz(:lorenz1)
lorenz2 = lorenz(:lorenz2)

connected = ODESystem([s ~ a + lorenz1.x
                       lorenz2.y ~ s
                       lorenz1.F ~ lorenz2.u
                       lorenz2.F ~ lorenz1.u],t,systems=[lorenz1,lorenz2])

# Reduced Flattened System

reduced_system = alias_elimination(connected; conservative=false)

@test setdiff(states(reduced_system), [
        a
        lorenz1.x
        lorenz1.y
        lorenz1.z
        lorenz2.x
        lorenz2.y
        lorenz2.z
       ]) |> isempty

@test setdiff(parameters(reduced_system), [
        lorenz1.σ
        lorenz1.ρ
        lorenz1.β
        lorenz2.σ
        lorenz2.ρ
        lorenz2.β
       ]) |> isempty

reduced_eqs = [
               0 ~ a + lorenz1.x - (lorenz2.y)
               D(lorenz1.x) ~ lorenz2.x + lorenz2.y + lorenz1.σ*((lorenz1.y) - (lorenz1.x)) - (lorenz2.z)
               D(lorenz1.y) ~ lorenz1.x*(lorenz1.ρ - (lorenz1.z)) - ((lorenz1.x) + (lorenz1.y) - (lorenz1.z))
               D(lorenz1.z) ~ lorenz1.x*lorenz1.y - (lorenz1.β*(lorenz1.z))
               D(lorenz2.x) ~ lorenz1.x + lorenz1.y + lorenz2.σ*((lorenz2.y) - (lorenz2.x)) - (lorenz1.z)
               D(lorenz2.y) ~ lorenz2.x*(lorenz2.ρ - (lorenz2.z)) - ((lorenz2.x) + (lorenz2.y) - (lorenz2.z))
               D(lorenz2.z) ~ lorenz2.x*lorenz2.y - (lorenz2.β*(lorenz2.z))
              ]

test_equal.(equations(reduced_system), reduced_eqs)

observed_eqs = [
                s ~ a + lorenz1.x
                lorenz1.u ~ lorenz1.x + lorenz1.y - lorenz1.z
                lorenz2.u ~ lorenz2.x + lorenz2.y - lorenz2.z
                lorenz2.F ~ lorenz1.u
                lorenz1.F ~ lorenz2.u
               ]
test_equal.(observed(reduced_system), observed_eqs)

pp = [
      lorenz1.σ => 10
      lorenz1.ρ => 28
      lorenz1.β => 8/3
      lorenz2.σ => 10
      lorenz2.ρ => 28
      lorenz2.β => 8/3
     ]
u0 = [
      a         => 1.0
      lorenz1.x => 1.0
      lorenz1.y => 0.0
      lorenz1.z => 0.0
      lorenz2.x => 1.0
      lorenz2.y => 0.0
      lorenz2.z => 0.0
     ]
prob1 = ODEProblem(reduced_system, u0, (0.0, 100.0), pp)
solve(prob1, Rodas5())

# issue #578

let
    @variables t x(t) y(t) z(t)
    D = Differential(t)
    eqs = [
        D(x) ~ x + y
        x ~ y
    ]
    sys = ODESystem(eqs, t)
    asys = alias_elimination(sys)

    test_equal.(asys.eqs, [D(x) ~ 2x])
    test_equal.(asys.observed, [y ~ x])
end

# issue #724 and #716
let
    @parameters t
    D = Differential(t)
    @variables x(t), u(t), y(t)
    @parameters a, b, c, d
    ol = ODESystem([D(x) ~ a * x + b * u; y ~ c * x + d * u], t, name=:ol)
    @variables u_c(t), y_c(t)
    @parameters k_P
    pc = ODESystem(Equation[u_c ~ k_P * y_c], t, name=:pc)
    connections = [
        ol.u ~ pc.u_c
        pc.y_c ~ ol.y
    ]
    connected = ODESystem(connections, t, systems=[ol, pc])
    @test equations(connected) isa Vector{Equation}
    reduced_sys = alias_elimination(connected)
    ref_eqs = [
               D(ol.x) ~ ol.a*ol.x + ol.b*pc.u_c
               0 ~ ol.c*ol.x + ol.d*pc.u_c - ol.y
               0 ~ pc.k_P*ol.y - pc.u_c
              ]
    @test ref_eqs == equations(reduced_sys)
end

# NonlinearSystem
@parameters t
@variables u1(t) u2(t) u3(t) u4(t) u5(t)
@parameters p
eqs = [
       2u1 ~ 3u5
       u2 ~ u1
       u3 ~ 2u1 - u2
       u4 ~ u2 + u3^p
       u5 ~ u4^2 - u1
      ]
sys = NonlinearSystem(eqs, [u1, u2, u3, u4, u5], [p])
reducedsys = alias_elimination(sys)
@test observed(reducedsys) == [u1 ~ 3/2 * u5]

u0 = [
      u1 => 1
      u2 => 1
      u3 => 0.3
      u4 => 0.6
      u5 => 2/3
     ]
pp = [2]
nlprob = NonlinearProblem(reducedsys, u0, pp)
reducedsol = solve(nlprob, NewtonRaphson())
residual = fill(100.0, 4)
nlprob.f(residual, reducedsol.u, pp)
@test all(x->abs(x) < 1e-5, residual)
