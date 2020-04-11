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

connected1 = ODESystem(connnectedeqs,t,[a],[α],systems=[lorenz1,lorenz2],name=:connected1)

@variables lorenz1′x(t) lorenz1′y(t) lorenz1′z(t) lorenz2′x(t) lorenz2′y(t) lorenz2′z(t)
@parameters lorenz1′σ lorenz1′ρ lorenz1′β lorenz2′σ lorenz2′ρ lorenz2′β

eqs_flat = [D(a) ~ a*lorenz1′x,
            D(lorenz1′x) ~ lorenz1′σ*(lorenz1′y-lorenz1′x),
            D(lorenz1′y) ~ lorenz1′x*(lorenz1′ρ-lorenz1′z)-lorenz1′y,
            0 ~ lorenz1′x + lorenz1′y + lorenz1′β*lorenz1′z,
            D(lorenz2′x) ~ lorenz2′σ*(lorenz2′y-lorenz2′x),
            D(lorenz2′y) ~ lorenz2′x*(lorenz2′ρ-lorenz2′z)-lorenz2′y,
            0 ~ lorenz2′x + lorenz2′y + lorenz2′β*lorenz2′z]

@test [x.name for x in states(connected1)] == [:a,:lorenz1′x,:lorenz1′y,:lorenz1′z,:lorenz2′x,:lorenz2′y,:lorenz2′z]
@test [x.name for x in parameters(connected1)] == [:α,:lorenz1′σ,:lorenz1′ρ,:lorenz1′β,:lorenz2′σ,:lorenz2′ρ,:lorenz2′β]
@test eqs_flat == equations(connected1)

connected2 = ODESystem(connnectedeqs,t,[a],[α],systems=[lorenz1,lorenz2],name=:connected2)

@parameters γ
@variables g(t)
connnectedeqs2 = [D(g) ~ g*states(connected1,lorenz1,:x)]
doublelevel = ODESystem(connnectedeqs2,t,[g],[γ],systems=[connected1,connected2],name=:doublelevel)

@test [x.name for x in states(doublelevel)] == [:g,
                                                :connected1′a,:connected1′lorenz1′x,:connected1′lorenz1′y,:connected1′lorenz1′z,:connected1′lorenz2′x,:connected1′lorenz2′y,:connected1′lorenz2′z,
                                                :connected2′a,:connected2′lorenz1′x,:connected2′lorenz1′y,:connected2′lorenz1′z,:connected2′lorenz2′x,:connected2′lorenz2′y,:connected2′lorenz2′z]
@test [x.name for x in parameters(doublelevel)] == [:γ,
                                                   :connected1′α,:connected1′lorenz1′σ,:connected1′lorenz1′ρ,:connected1′lorenz1′β,:connected1′lorenz2′σ,:connected1′lorenz2′ρ,:connected1′lorenz2′β,
                                                   :connected2′α,:connected2′lorenz1′σ,:connected2′lorenz1′ρ,:connected2′lorenz1′β,:connected2′lorenz2′σ,:connected2′lorenz2′ρ,:connected2′lorenz2′β]

@variables connected1′a(t) connected1′lorenz1′x(t) connected1′lorenz1′y(t) connected1′lorenz1′z(t) connected1′lorenz2′x(t) connected1′lorenz2′y(t) connected1′lorenz2′z(t)
@variables connected2′a(t) connected2′lorenz1′x(t) connected2′lorenz1′y(t) connected2′lorenz1′z(t) connected2′lorenz2′x(t) connected2′lorenz2′y(t) connected2′lorenz2′z(t)
@parameters connected1′α connected1′lorenz1′σ connected1′lorenz1′ρ connected1′lorenz1′β connected1′lorenz2′σ connected1′lorenz2′ρ connected1′lorenz2′β
@parameters connected2′α connected2′lorenz1′σ connected2′lorenz1′ρ connected2′lorenz1′β connected2′lorenz2′σ connected2′lorenz2′ρ connected2′lorenz2′β

eqs_flat2 = [D(g) ~ g*connected1′lorenz1′x,
            D(connected1′a) ~ connected1′a*connected1′lorenz1′x,
            D(connected1′lorenz1′x) ~ connected1′lorenz1′σ*(connected1′lorenz1′y-connected1′lorenz1′x),
            D(connected1′lorenz1′y) ~ connected1′lorenz1′x*(connected1′lorenz1′ρ-connected1′lorenz1′z)-connected1′lorenz1′y,
            0 ~ connected1′lorenz1′x + connected1′lorenz1′y + connected1′lorenz1′β*connected1′lorenz1′z,
            D(connected1′lorenz2′x) ~ connected1′lorenz2′σ*(connected1′lorenz2′y-connected1′lorenz2′x),
            D(connected1′lorenz2′y) ~ connected1′lorenz2′x*(connected1′lorenz2′ρ-connected1′lorenz2′z)-connected1′lorenz2′y,
            0 ~ connected1′lorenz2′x + connected1′lorenz2′y + connected1′lorenz2′β*connected1′lorenz2′z,
            D(connected2′a) ~ connected2′a*connected2′lorenz1′x,
            D(connected2′lorenz1′x) ~ connected2′lorenz1′σ*(connected2′lorenz1′y-connected2′lorenz1′x),
            D(connected2′lorenz1′y) ~ connected2′lorenz1′x*(connected2′lorenz1′ρ-connected2′lorenz1′z)-connected2′lorenz1′y,
            0 ~ connected2′lorenz1′x + connected2′lorenz1′y + connected2′lorenz1′β*connected2′lorenz1′z,
            D(connected2′lorenz2′x) ~ connected2′lorenz2′σ*(connected2′lorenz2′y-connected2′lorenz2′x),
            D(connected2′lorenz2′y) ~ connected2′lorenz2′x*(connected2′lorenz2′ρ-connected2′lorenz2′z)-connected2′lorenz2′y,
            0 ~ connected2′lorenz2′x + connected2′lorenz2′y + connected2′lorenz2′β*connected2′lorenz2′z]

@test eqs_flat2 == equations(doublelevel)

M = Array(I,15,15)
M[5,5] = false
M[8,8] = false
M[12,12] = false
M[15,15] = false
@test calculate_massmatrix(doublelevel) == M

jac = [connected1′lorenz1′x 0 g zeros(1,12)
       zeros(7,1) calculate_jacobian(connected1) zeros(Expression,7,7)
       zeros(Expression,7,8) calculate_jacobian(connected2)]

jac2 = [connected1′lorenz1′x 0 g zeros(1,12)
        zeros(7,1) ModelingToolkit.namespace_operation.(calculate_jacobian(connected1),connected1.name,:t) zeros(Expression,7,7)
        zeros(Expression,7,8) ModelingToolkit.namespace_operation.(calculate_jacobian(connected2),connected2.name,:t)]

@test all(isequal.(calculate_jacobian(doublelevel),jac2))
