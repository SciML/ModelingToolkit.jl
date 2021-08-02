using ModelingToolkit, Unitful
using Test
MT = ModelingToolkit
@parameters τ [unit=u"ms"]
@variables t [unit=u"ms"] E(t) [unit=u"kJ"] P(t) [unit=u"MW"]
D = Differential(t)

@test MT.vartype(t) == u"ms"
@test MT.vartype(E) == u"kJ"
@test MT.vartype(τ) == u"ms"

eqs = [D(E) ~ P-E/τ ]
sys = ODESystem(eqs)

@test MT.instantiate(eqs[1].lhs) == 1.0u"MW"
@test MT.instantiate(eqs[1].rhs) == 1.0u"MW"
@test MT.validate(eqs[1])
@test MT.validate(sys)

@test MT.instantiate(0.5) == 1.0
@test MT.instantiate(t) == 1.0u"ms"
@test MT.instantiate(P) == 1.0u"MW"
@test MT.instantiate(τ) == 1.0u"ms"

@test MT.instantiate(τ^-1) == 1/u"ms"
@test MT.instantiate(D(E)) == 1.0u"MW"
@test MT.instantiate(E/τ) == 1.0u"MW" 
@test MT.instantiate(2*P) == 1.0u"MW"
@test MT.instantiate(t/τ) == 1.0
@test MT.instantiate(P-E/τ)/1.0u"MW" == 1.0

@test MT.instantiate(1.0^(t/τ)) == 1.0
@test MT.instantiate(exp(t/τ)) == 1.0
@test MT.instantiate(sin(t/τ)) == 1.0
@test MT.instantiate(sin(1u"rad")) == 1.0
@test MT.instantiate(t^2) == 1.0u"ms"^2

@test !MT.validate(E^1.5~ E^(t/τ))
@test MT.validate(E^(t/τ)~ E^(t/τ))

sys = ODESystem(eqs,t,[P,E],[τ])
@test MT.validate(sys)

@test !MT.validate(D(D(E))~P)
@test !MT.validate(0~P+E*τ)
@test_logs (:warn,) MT.validate(0 ~ P + E*τ)
@test_logs (:warn,) MT.validate(P + E*τ ~ 0)
@test_logs (:warn,) MT.validate(P ~ 0)

@variables x y z u
@parameters σ ρ β
eqs = [0 ~ σ*(y-x)]
@test MT.validate(eqs) #should cope with unit-free

@variables t x[1:3,1:3](t) #should cope with arrays
D = Differential(t)
eqs = D.(x) .~ x
ODESystem(eqs)