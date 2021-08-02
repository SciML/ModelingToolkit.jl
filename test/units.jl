using ModelingToolkit, Unitful, OrdinaryDiffEq
using Test
MT = ModelingToolkit
@parameters τ [unit = u"ms"]
@variables t [unit = u"ms"] E(t) [unit = u"kJ"] P(t) [unit = u"MW"]
D = Differential(t)

@test MT.vartype(t) == u"ms"
@test MT.vartype(E) == u"kJ"
@test MT.vartype(τ) == u"ms"


@test MT.instantiate(0.5) == 1.0
@test MT.instantiate(t) == 1.0u"ms"
@test MT.instantiate(P) == 1.0u"MW"
@test MT.instantiate(τ) == 1.0u"ms"

@test MT.instantiate(τ^-1) == 1/u"ms"
@test MT.instantiate(D(E)) == 1.0u"MW"
@test MT.instantiate(E/τ) == 1.0u"MW"
@test MT.instantiate(2*P) == 1.0u"MW"
@test MT.instantiate(t/τ) == 1.0
@test MT.instantiate(P - E/τ)/1.0u"MW" == 1.0

@test MT.instantiate(1.0^(t/τ)) == 1.0
@test MT.instantiate(exp(t/τ)) == 1.0
@test MT.instantiate(sin(t/τ)) == 1.0
@test MT.instantiate(sin(1u"rad")) == 1.0
@test MT.instantiate(t^2) == 1.0u"ms"^2

@test !MT.validate(E^1.5 ~ E^(t/τ))
@test MT.validate(E^(t/τ) ~ E^(t/τ))

eqs = [D(E) ~ P - E/τ
        0.0u"MW" ~ P]
@test MT.instantiate(eqs[1].lhs) == 1.0u"MW"
@test MT.instantiate(eqs[1].rhs) == 1.0u"MW"
@test MT.validate(eqs[1])
@test MT.validate(eqs[2])
sys = ODESystem(eqs)
sys = ODESystem(eqs, t, [P, E], [τ])
@test MT.validate(sys)

@test !MT.validate(D(D(E)) ~ P)
@test !MT.validate(0 ~ P + E*τ)
@test_logs (:warn,) MT.validate(0 ~ P + E*τ)
@test_logs (:warn,) MT.validate(P + E*τ ~ 0)
@test_logs (:warn,) MT.validate(P ~ 0)

#Unit-free
@variables x y z u
@parameters σ ρ β
eqs = [0 ~ σ*(y - x)]
@test MT.validate(eqs)

#Array variables
@variables t x[1:3,1:3](t)
D = Differential(t)
eqs = D.(x) .~ x
ODESystem(eqs)

# Array ops
using Symbolics: unwrap, wrap
using LinearAlgebra
@variables t
sts = @variables x[1:3](t) y(t)
ps = @parameters p[1:3] = [1, 2, 3]
D = Differential(t)
eqs = [
       collect(D.(x) ~ x)
       D(y) ~ norm(x)*y
      ]
ODESystem(eqs, t, [sts...;], [ps...;])

#= Not supported yet b/c iterate doesn't work on unitful array
# Array ops with units
@variables t [unit =u"s"]
sts = @variables x[1:3](t) [unit = u"kg"] y(t) [unit = u"kg"]
ps = @parameters b [unit = u"s"^-1]
D = Differential(t)
eqs = [
       collect(D.(x) ~ b*x)
       D(y) ~ b*norm(x)
      ]
ODESystem(eqs, t, [sts...;], [ps...;])

#Array variables with units
@variables t [unit = u"s"] x[1:3,1:3](t) [unit = u"kg"] 
@parameters a [unit = u"s"^-1]
D = Differential(t)
eqs = D.(x) .~ a*x
ODESystem(eqs)
=#

#Difference equation with units
@parameters t [unit = u"s"] a [unit = u"s"^-1]
@variables x(t) [unit = u"kg"]
δ = Differential(t)
D = Difference(t; dt = 0.1u"s")
eqs = [
    δ(x) ~ a*x 
]
de = ODESystem(eqs, t, [x, y], [a])


@parameters t
@variables y[1:3](t)
@parameters k[1:3]
D = Differential(t)

eqs = [D(y[1]) ~ -k[1]*y[1] + k[3]*y[2]*y[3],
       D(y[2]) ~  k[1]*y[1] - k[3]*y[2]*y[3] - k[2]*y[2]^2,
       0 ~  y[1] + y[2] + y[3] - 1]

sys = ODESystem(eqs,t,y,k)