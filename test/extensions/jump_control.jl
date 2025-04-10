using ModelingToolkit
using JuMP, InfiniteOpt
using DiffEqDevTools, DiffEqBase
using SimpleDiffEq
using OrdinaryDiffEqSDIRK
using Ipopt
using BenchmarkTools
const M = ModelingToolkit

@testset "ODE Solution, no cost" begin
    # Test solving without anything attached.
    @parameters α=1.5 β=1.0 γ=3.0 δ=1.0
    M.@variables x(..) y(..)
    t = M.t_nounits
    D = M.D_nounits

    eqs = [D(x(t)) ~ α * x(t) - β * x(t) * y(t),
        D(y(t)) ~ -γ * y(t) + δ * x(t) * y(t)]

    tspan = (0.0, 1.0)
    u0map = [x(t) => 4.0, y(t) => 2.0]
    parammap = [α => 1.5, β => 1.0, γ => 3.0, δ => 1.0]
    @mtkbuild sys = ODESystem(eqs, t)

    # Test explicit method.
    jprob = JuMPControlProblem(sys, u0map, tspan, parammap, dt = 0.01)
    @test num_constraints(jprob.model) == 2 # initials
    jsol = solve(jprob, Ipopt.Optimizer, :RK4)
    oprob = ODEProblem(sys, u0map, tspan, parammap)
    osol = solve(oprob, SimpleRK4(), dt = 0.01)
    @test jsol.sol.u ≈ osol.u

    # Implicit method.
    jsol2 = @btime solve($jprob, Ipopt.Optimizer, :ImplicitEuler, silent = true) # 63.031 ms, 26.49 MiB
    osol2 = @btime solve($oprob, ImplicitEuler(), dt = 0.01, adaptive = false) # 129.375 μs, 61.91 KiB
    @test ≈(jsol2.sol.u, osol2.u, rtol = 0.001)
    iprob = InfiniteOptControlProblem(sys, u0map, tspan, parammap, dt = 0.01)
    isol = @btime solve(
        $iprob, Ipopt.Optimizer, derivative_method = FiniteDifference(Backward()), silent = true) # 11.540 ms, 4.00 MiB

    # With a constraint
    u0map = Pair[]
    guess = [x(t) => 4.0, y(t) => 2.0]
    constr = [x(0.6) ~ 3.5, x(0.3) ~ 7.0]
    @mtkbuild lksys = ODESystem(eqs, t; constraints = constr)

    jprob = JuMPControlProblem(lksys, u0map, tspan, parammap; guesses = guess, dt = 0.01)
    @test num_constraints(jprob.model) == 2
    jsol = @btime solve($jprob, Ipopt.Optimizer, :Tsitouras5, silent = true) # 12.190 s, 9.68 GiB
    sol = jsol.sol
    @test sol(0.6)[1] ≈ 3.5
    @test sol(0.3)[1] ≈ 7.0

    iprob = InfiniteOptControlProblem(
        lksys, u0map, tspan, parammap; guesses = guess, dt = 0.01)
    isol = @btime solve(
        $iprob, Ipopt.Optimizer, derivative_method = OrthogonalCollocation(3), silent = true) # 48.564 ms, 9.58 MiB
    sol = isol.sol
    @test sol(0.6)[1] ≈ 3.5
    @test sol(0.3)[1] ≈ 7.0
end

#@testset "Optimal control: bees" begin
#    # Example from Lawrence Evans' notes
#    M.@variables w(..) q(..) 
#    M.@parameters α(t) [bounds = [0, 1]] b c μ s ν
#
#    tspan = (0, 4)
#    eqs = [D(w(t)) ~ -μ*w(t) + b*s*α*w(t),
#           D(q(t)) ~ -ν*q(t) + c*(1 - α)*s*w(t)]
#    costs = [-q(tspan[2])]
#   
#    @mtkbuild beesys = ODESystem(eqs, t; costs)
#    u0map = [w(0) => 40, q(0) => 2]
#    pmap = [b => 1, c => 1, μ => 1, s => 1, ν => 1]
#
#    jprob = JuMPControlProblem(beesys, u0map, tspan, pmap)
#    jsol = solve(jprob, Ipopt.Optimizer, :Tsitouras5)
#    control_sol = jsol.control_sol
#    # Bang-bang control
#end
#
#@testset "Constrained optimal control problems" begin
#end
