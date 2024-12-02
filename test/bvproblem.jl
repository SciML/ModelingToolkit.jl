using BoundaryValueDiffEq, OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters σ = 10. ρ = 28 β = 8/3
@variables x(t) = 1 y(t) = 0 z(t) = 0

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]

u0map = [:x => 1., :y => 0., :z => 0.]
parammap = [:ρ => 28., :β => 8/3, :σ => 10.]
tspan = (0., 10.)

@mtkbuild lorenz = ODESystem(eqs, t)

bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(lorenz, u0map, tspan, parammap)
sol = solve(bvp, MIRK4(), dt = 0.1);

bvp2 = SciMLBase.BVProblem{false, SciMLBase.AutoSpecialize}(lorenz, u0map, tspan, parammap)
sol2 = solve(bvp, MIRK4(), dt = 0.1);

op = ODEProblem(lorenz, u0map, tspan, parammap)
osol = solve(op)

@test sol.u[end] ≈ osol.u[end]
@test sol2.u[end] ≈ osol.u[end]
@test sol.u[1] == [1., 0., 0.]
@test sol2.u[1] == [1., 0., 0.]

### Testing on pendulum

@parameters g = 9.81 L = 1. 
@variables θ(t) = π/2 

eqs = [D(D(θ)) ~ -(g / L) * sin(θ)]

@mtkbuild pend = ODESystem(eqs, t)

u0map = [θ => π/2, D(θ) => π/2]
parammap = [:L => 2.]
tspan = (0., 10.)

bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap)
sol = solve(bvp, MIRK4(), dt = 0.05);

bvp2 = SciMLBase.BVProblem{false, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap)
sol2 = solve(bvp2, MIRK4(), dt = 0.05);

osol = solve(pend)

@test sol.u[end] ≈ osol.u[end]
@test sol.u[1] == [π/2, π/2]
@test sol2.u[end] ≈ osol.u[end]
@test sol2.u[1] == [π/2, π/2]
