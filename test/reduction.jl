using ModelingToolkit, OrdinaryDiffEq, Test

@parameters t σ ρ β
@variables x(t) y(t) z(t) a(t)
@derivatives D'~t

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ a*y - β*z,
       0 ~ x - a]

lorenz1 = ODESystem(eqs,t,[x,y,z,a],[σ,ρ,β],name=:lorenz1)

lorenz1_aliased = alias_elimination(lorenz1)
length(equations(lorenz1_reduced)) = 3
length(states(lorenz1_reduced)) = 3

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]
test_lorenz1_aliased = ODESystem(eqs,t,[x,y,z],[σ,ρ,β],observed=[a ~ x],name=:lorenz1)

# Multi-System Reduction

eqs1 = [D(x) ~ σ*(y-x) + F,
       D(y) ~ x*(ρ-z)-u,
       D(z) ~ x*y - β*z]
aliases = [u ~ x + y - z]
lorenz1 = ODESystem(eqs1,pins=[F],observed=aliases,name=:lorenz1)

eqs2 = [D(x) ~ F,
       D(y) ~ x*(ρ-z)-x,
       D(z) ~ x*y - β*z]
aliases2 = [u ~ x - y - z]
lorenz2 = ODESystem(eqs2,pins=[F],observed=aliases2,name=:lorenz2)

connections = [lorenz1.F ~ lorenz2.u,
               lorenz2.F ~ lorenz1.u]
connected = ODESystem([0 ~ a + lorenz1.x - lorenz2.y],t,[a],[],observed=connections,systems=[lorenz1,lorenz2])

# Reduced Unflattened System

connections2 = [lorenz1.F ~ lorenz2.u,
                lorenz2.F ~ lorenz1.u,
                a ~ -lorenz1.x + lorenz2.y]
connected = ODESystem(Equation[],t,[],[],observed=connections2,systems=[lorenz1,lorenz2])

# Reduced Flattened System

flattened_system = flatten(connected)
flatten(sys::AbstractSystem) = ODESystem(equations(sys),states(sys),parameters(sys),independent_variable(sys))

aliased_flattened_system = alias_elimination(flattened_system)

states(reduced_flattened_system) == [
        lorenz1.x
        lorenz1.y
        lorenz1.z
        lorenz2.x
        lorenz2.y
        lorenz2.z
]

parameters(reduced_flattened_system) == [
        lorenz1.σ
        lorenz1.ρ
        lorenz1.β
        lorenz2.σ
        lorenz2.ρ
        lorenz2.β
]

equations(reduced_flattened_system) == [
       D(lorenz1.x) ~ lorenz1.σ*(lorenz1.y-lorenz1.x) + lorenz2.x - lorenz2.y - lorenz2.z,
       D(lorenz1.y) ~ lorenz1.x*(ρ-z)-lorenz1.x - lorenz1.y + lorenz1.z,
       D(lorenz1.z) ~ lorenz1.x*lorenz1.y - lorenz1.β*lorenz1.z
       D(lorenz2.x) ~ lorenz1.x + lorenz1.y - lorenz1.z,
       D(lorenz2.y) ~ lorenz2.x*(lorenz2.ρ-lorenz2.z)-lorenz2.x,
       D(lorenz2.z) ~ lorenz2.x*lorenz2.y - lorenz2.β*lorenz2.z]

observed(reduced_flattened_system) == [
        lorenz1.F ~ lorenz2.u
        lorenz2.F ~ lorenz1.u
        a ~ -lorenz1.x + lorenz2.y
]
