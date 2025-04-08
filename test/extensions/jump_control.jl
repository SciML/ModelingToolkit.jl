using ModelingToolkit
using JuMP, InfiniteOpt
using DiffEqDevTools, DiffEqBase
using SimpleDiffEq
using HiGHS
const M = ModelingToolkit

@testset "ODE Solution, no cost" begin
    # Test solving without anything attached.
    @parameters α=1.5 β=1.0 γ=3.0 δ=1.0
    M.@variables x(..) y(..)
    t = M.t_nounits; D = M.D_nounits

    eqs = [D(x(t)) ~ α * x(t) - β * x(t) * y(t),
        D(y(t)) ~ -γ * y(t) + δ * x(t) * y(t)]

    tspan = (0.0, 1.0)
    u0map = [x(t) => 4.0, y(t) => 2.0]
    parammap = [α => 7.5, β => 4, γ => 8.0, δ => 5.0]
    @mtkbuild sys = ODESystem(eqs, t)

    # Test explicit method.
    jprob = JuMPControlProblem(sys, u0map, tspan, parammap, dt = 0.01)
    @test num_constraints(jprob.model) == 2 # initials
    jsol = solve(jprob, Ipopt.Optimizer, :Tsitouras5)
    oprob = ODEProblem(sys, u0map, tspan, parammap)
    osol = solve(oprob, SimpleTsit5(), adaptive = false)
    @test jsol.sol.u ≈ osol.u

    # Implicit method.
    jsol2 = solve(prob, Ipopt.Optimizer, :RK4)
    osol2 = solve(oprob, SimpleRK4(), adaptive = false)
    @test jsol2.sol.u ≈ osol2.u

    # With a constraint
    @parameters α=1.5 β=1.0 γ=3.0 δ=1.0
    @variables x(..) y(..)

    eqs = [D(x(t)) ~ α * x(t) - β * x(t) * y(t),
        D(y(t)) ~ -γ * y(t) + δ * x(t) * y(t)]

    u0map = []
    tspan = (0.0, 1.0)
    guess = [x(t) => 4.0, y(t) => 2.0]
    constr = [x(0.6) ~ 3.5, x(0.3) ~ 7.0]
    @mtkbuild lksys = ODESystem(eqs, t; constraints = constr)

    jprob = JuMPControlProblem(sys, u0map, tspan, parammap; guesses, dt = 0.01)
    @test num_constraints(jprob.model) == 2 == num_variables(jprob.model) == 2
    jsol = solve(prob, HiGHS.Optimizer, :Tsitouras5)
    sol = jsol.sol
    @test sol(0.6)[1] ≈ 3.5
    @test sol(0.3)[1] ≈ 7.0
end

@testset "Optimal control problem" begin
end
