using BoundaryValueDiffEq, OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

solvers = [MIRK4, RadauIIa5, LobattoIIIa3]

@parameters α=7.5 β=4.0 γ=8.0 δ=5.0
@variables x(t)=1.0 y(t)=2.0

eqs = [D(x) ~ α * x - β * x * y,
    D(y) ~ -γ * y + δ * x * y]

u0map = [:x => 1.0, :y => 2.0]
parammap = [:α => 7.5, :β => 4, :γ => 8.0, :δ => 5.0]
tspan = (0.0, 10.0)

@mtkbuild lotkavolterra = ODESystem(eqs, t)
op = ODEProblem(lotkavolterra, u0map, tspan, parammap)
osol = solve(op, Vern9())

bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(
    lotkavolterra, u0map, tspan, parammap; eval_expression = true)

for solver in solvers
    sol = solve(bvp, solver(), dt = 0.01)
    @test isapprox(sol.u[end], osol.u[end]; atol = 0.01)
    @test sol.u[1] == [1.0, 2.0]
end

# Test out of place
bvp2 = SciMLBase.BVProblem{false, SciMLBase.AutoSpecialize}(
    lotkavolterra, u0map, tspan, parammap; eval_expression = true)

for solver in solvers
    sol = solve(bvp2, solver(), dt = 0.01)
    @test isapprox(sol.u[end], osol.u[end]; atol = 0.01)
    @test sol.u[1] == [1.0, 2.0]
end

### Testing on pendulum

@parameters g=9.81 L=1.0
@variables θ(t) = π / 2

eqs = [D(D(θ)) ~ -(g / L) * sin(θ)]

@mtkbuild pend = ODESystem(eqs, t)

u0map = [θ => π / 2, D(θ) => π / 2]
parammap = [:L => 1.0, :g => 9.81]
tspan = (0.0, 6.0)

op = ODEProblem(pend, u0map, tspan, parammap)
osol = solve(op, Vern9())

bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap)
for solver in solvers
    sol = solve(bvp, solver(), dt = 0.01)
    @test isapprox(sol.u[end], osol.u[end]; atol = 0.01)
    @test sol.u[1] == [π / 2, π / 2]
end

# Test out-of-place
bvp2 = SciMLBase.BVProblem{false, SciMLBase.FullSpecialize}(pend, u0map, tspan, parammap)

for solver in solvers
    sol = solve(bvp2, solver(), dt = 0.01)
    @test isapprox(sol.u[end], osol.u[end]; atol = 0.01)
    @test sol.u[1] == [π / 2, π / 2]
end
