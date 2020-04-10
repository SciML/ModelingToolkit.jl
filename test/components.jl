using ModelingToolkit, StaticArrays, LinearAlgebra
using DiffEqBase
using Test

# Define some variables
@parameters t σ ρ β
@variables x(t) y(t) z(t)
@derivatives D'~t

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       0 ~ x + y + z]

de1 = ODESystem(eqs,t,[x,y,z],[σ,ρ,β],name=:lorenz1)
de2 = ODESystem(eqs,t,[x,y,z],[σ,ρ,β],name=:lorenz2)

@parameters α
@variables a(t)
for

convert(DAESystem,sys)
convert(ODESystem,sys)

connected = ODESystems(connnectedeqs,[a],[α],systems=[de1,de2],name=:connectedlorenz)

deleteat!(de2.eqs,3)
de2
push!(de2.eqs,convert(ModelingToolkit.ODEExpr,D(z) ~ z^2)[2])

check_consistency(de2)

b = a*z
           [states(de1,:x) ~ α*states(de2,:z),
            D(a) ~ b*states(de1,:x)]

states(states) == [lorenz1′x,lorenz1′y,lorenz1′z,lorenz2′x,lorenz2′y,lorenz2′z,a]
parameters(connected) == [lorenz1′σ,lorenz1′ρ,lorenz1′β,lorenz2′σ,lorenz2′ρ,lorenz2′β,α]
