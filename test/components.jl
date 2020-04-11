using ModelingToolkit, StaticArrays, LinearAlgebra
using DiffEqBase
using Test

# Define some variables
@parameters t σ ρ β
@variables x(t) y(t) z(t)
@derivatives D'~t

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       0 ~ x + y + β*z]

lorenz1 = ODESystem(eqs,t,[x,y,z],[σ,ρ,β],name=:lorenz1)
lorenz2 = ODESystem(eqs,t,[x,y,z],[σ,ρ,β],name=:lorenz2)

@parameters α
@variables a(t)
connnectedeqs = [D(a) ~ a*states(lorenz1,:x)]

connected = ODESystem(connnectedeqs,t,[a],[α],systems=[lorenz1,lorenz2],name=:connectedlorenz)

@variables lorenz1′x(t) lorenz1′y(t) lorenz1′z(t) lorenz2′x(t) lorenz2′y(t) lorenz2′z(t)
@parameters lorenz1′σ lorenz1′ρ lorenz1′β lorenz2′σ lorenz2′ρ lorenz2′β

eqs_flat = [D(a) ~ a*lorenz1′x,
            D(lorenz1′x) ~ lorenz1′σ*(lorenz1′y-lorenz1′x),
            D(lorenz1′y) ~ lorenz1′x*(lorenz1′ρ-lorenz1′z)-lorenz1′y,
            0 ~ lorenz1′x + lorenz1′y + lorenz1′β*lorenz1′z,
            D(lorenz2′x) ~ lorenz2′σ*(lorenz2′y-lorenz2′x),
            D(lorenz2′y) ~ lorenz2′x*(lorenz2′ρ-lorenz2′z)-lorenz2′y,
            0 ~ lorenz2′x + lorenz2′y + lorenz2′β*lorenz2′z]

@test [x.name for x in states(connected)] == [:a,:lorenz1′x,:lorenz1′y,:lorenz1′z,:lorenz2′x,:lorenz2′y,:lorenz2′z]
@test [x.name for x in parameters(connected)] == [:α,:lorenz1′σ,:lorenz1′ρ,:lorenz1′β,:lorenz2′σ,:lorenz2′ρ,:lorenz2′β]
@test eqs_flat == equations(connected)
