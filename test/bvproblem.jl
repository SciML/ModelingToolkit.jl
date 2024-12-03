using BoundaryValueDiffEq, OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters α = 7.5 β = 4. γ = 8. δ = 5. 
@variables x(t) = 1. y(t) = 2. 

eqs = [D(x) ~ α*x - β*x*y,
       D(y) ~ -γ*y + δ*x*y]

u0map = [:x => 1., :y => 2.]
parammap = [:α => 7.5, :β => 4, :γ => 8., :δ => 5.]
tspan = (0., 10.)

@mtkbuild lotkavolterra = ODESystem(eqs, t)

bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(lotkavolterra, u0map, tspan, parammap)
sol = solve(bvp, MIRK4(), dt = 0.01);

bvp2 = SciMLBase.BVProblem{false, SciMLBase.AutoSpecialize}(lotkavolterra, u0map, tspan, parammap)
sol2 = solve(bvp, MIRK4(), dt = 0.01);

op = ODEProblem(lotkavolterra, u0map, tspan, parammap)
osol = solve(op, Vern9())

@test isapprox(sol.u[end],osol.u[end]; atol = 0.001)
@test isapprox(sol2.u[end],osol.u[end]; atol = 0.001)
@test sol.u[1] == [1., 2.]
@test sol2.u[1] == [1., 2.]

### Testing on pendulum

@parameters g = 9.81 L = 1. 
@variables θ(t) = π/2 

eqs = [D(D(θ)) ~ -(g / L) * sin(θ)]

@mtkbuild pend = ODESystem(eqs, t)

u0map = [θ => π/2, D(θ) => π/2]
parammap = [:L => 2., :g => 9.81]
tspan = (0., 10.)

bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap)
sol = solve(bvp, MIRK4(), dt = 0.01);

bvp2 = SciMLBase.BVProblem{false, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap)
sol2 = solve(bvp2, MIRK4(), dt = 0.01);

op = ODEProblem(pend, u0map, tspan, parammap)
osol = solve(op, Vern9())

@test sol.u[end] ≈ osol.u[end]
@test sol.u[1] == [π/2, π/2]
@test sol2.u[end] ≈ osol.u[end]
@test sol2.u[1] == [π/2, π/2]
