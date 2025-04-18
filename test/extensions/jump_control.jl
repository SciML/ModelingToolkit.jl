using ModelingToolkit
import JuMP, InfiniteOpt
using DiffEqDevTools, DiffEqBase
using SimpleDiffEq
using OrdinaryDiffEqSDIRK
using Ipopt
using BenchmarkTools
using CairoMakie
const M = ModelingToolkit

@testset "ODE Solution, no cost" begin
    # Test solving without anything attached.
    @parameters α=1.5 β=1.0 γ=3.0 δ=1.0
    @variables x(..) y(..)
    t = M.t_nounits
    D = M.D_nounits

    eqs = [D(x(t)) ~ α * x(t) - β * x(t) * y(t),
        D(y(t)) ~ -γ * y(t) + δ * x(t) * y(t)]

    @mtkbuild sys = ODESystem(eqs, t)
    tspan = (0.0, 1.0)
    u0map = [x(t) => 4.0, y(t) => 2.0]
    parammap = [α => 1.5, β => 1.0, γ => 3.0, δ => 1.0]

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

    # Test whole-interval constraints
    constr = [x(t) > 3, y(t) > 4]
    @mtkbuild lksys = ODESystem(eqs, t; constraints = constr)
    iprob = InfiniteOptControlProblem(
        lksys, u0map, tspan, parammap; guesses = guess, dt = 0.01)
    isol = @btime solve(
        $iprob, Ipopt.Optimizer, derivative_method = OrthogonalCollocation(3), silent = true) # 48.564 ms, 9.58 MiB
    sol = isol.sol
    @test all(u -> u .> [3, 4], sol.u)

    jprob = JuMPControlProblem(lksys, u0map, tspan, parammap; guesses = guess, dt = 0.01)
    jsol = @btime solve($jprob, Ipopt.Optimizer, :RadauIA3, silent = true) # 12.190 s, 9.68 GiB
    sol = jsol.sol
    @test all(u -> u .> [3, 4], sol.u)
end

@testset "Linear systems" begin
    function is_bangbang(input_sol, lbounds, ubounds, rtol = 1e-4)
        bangbang = true
        for v in 1:length(input_sol.u[1]) - 1
            all(i -> ≈(i[v], bounds[v]; rtol) || ≈(i[v], bounds[u]; rtol), input_sol.u) || (bangbang = false)
        end
        bangbang
    end

    # Double integrator
    t = M.t_nounits
    D = M.D_nounits
    @variables x(..) [bounds = (0., 0.25)] v(..)
    @variables u(t) [bounds = (-1., 1.), input = true]
    constr = [v(1.0) ~ 0.0]
    cost = [-x(1.0)] # Maximize the final distance.
    @named block = ODESystem([D(x(t)) ~ v(t), D(v(t)) ~ u], t; costs = cost, constraints = constr)
    block, input_idxs = structural_simplify(block, ([u],[]))

    u0map = [x(t) => 0., v(t) => 0.]
    tspan = (0., 1.)
    parammap = [u => 0.]
    jprob = JuMPControlProblem(block, u0map, tspan, parammap; dt = 0.01)
    jsol = solve(jprob, Ipopt.Optimizer, :Verner8)
    # Linear systems have bang-bang controls
    @test is_bangbang(jsol.input_sol, [-1.], [1.])
    # Test reached final position.
    @test ≈(jsol.sol.u[end][1], 0.25, rtol = 1e-5)

    iprob = InfiniteOptControlProblem(block, u0map, tspan, parammap; dt = 0.01)
    isol = solve(iprob, Ipopt.Optimizer; silent = true)
    @test is_bangbang(isol.input_sol, [-1.], [1.])
    @test ≈(isol.sol.u[end][1], 0.25, rtol = 1e-5)

    ###################
    ### Bee example ###
    ###################
    # From Lawrence Evans' notes
    @variables w(..) q(..) α(t) [input = true, bounds = (0, 1)]
    @parameters b c μ s ν

    tspan = (0, 4)
    eqs = [D(w(t)) ~ -μ*w(t) + b*s*α*w(t),
           D(q(t)) ~ -ν*q(t) + c*(1 - α)*s*w(t)]
    costs = [-q(tspan[2])]
   
    @named beesys = ODESystem(eqs, t; costs)
    beesys, input_idxs = structural_simplify(beesys, ([α],[]))
    u0map = [w(t) => 40, q(t) => 2]
    pmap = [b => 1, c => 1, μ => 1, s => 1, ν => 1, α => 1]

    jprob = JuMPControlProblem(beesys, u0map, tspan, pmap, dt = 0.01)
    jsol = solve(jprob, Ipopt.Optimizer, :Tsitouras5)
    @test is_bangbang(jsol.input_sol, [0.], [1.])
    iprob = InfiniteOptControlProblem(beesys, u0map, tspan, pmap, dt = 0.01)
    isol = solve(jprob, Ipopt.Optimizer, :Tsitouras5)
    @test is_bangbang(isol.input_sol, [0.], [1.])
end

@testset "Rocket launch" begin
    t = M.t_nounits
    D = M.D_nounits

    @parameters h_c m₀ h₀ g₀ D_c c Tₘ m_c
    @variables h(..) v(..) m(..) [bounds = (m_c, 1)] T(..) [input = true, bounds = (0, Tₘ)]
    drag(h, v) = D_c * v^2 * exp(-h_c * (h - h₀) / h₀)
    gravity(h) = g₀ * (h₀ / h)

    eqs = [D(h(t)) ~ v(t), 
           D(v(t)) ~ (T(t) - drag(h(t), v(t))) / m(t) - gravity(h(t)),
           D(m(t)) ~ -T(t) / c]

    (ts, te) = (0., 0.2)
    costs = [-h(te)]
    constraints = [T(te) ~ 0]
    @named rocket = ODESystem(eqs, t; costs, constraints)
    rocket, input_idxs = structural_simplify(rocket, ([T(t)], []))

    u0map = [h(t) => h₀, m(t) => m₀, v(t) => 0]
    pmap = [g₀ => 1, m₀ => 1.0, h_c => 500, c => 0.5*√(g₀*h₀), D_c => 0.5 * 620 * m₀/g₀, Tₘ => 3.5*g₀*m₀, T(t) => 0., h₀ => 1, m_c => 0.6]
    jprob = JuMPControlProblem(rocket, u0map, (ts, te), pmap; dt = 0.005, cse = false)
    jsol = solve(jprob, Ipopt.Optimizer, :RadauIA3)
    @test jsol.sol.u[end][1] > 1.012
end

@testset "Free final time problem" begin
    t = M.t_nounits
    D = M.D_nounits

    @variables x(..) u(..) [input = true, bounds = (0,1)]
    @parameters tf
    eqs = [D(x(t)) ~ -2 + 0.5*u(t)]
    # Integral cost function
    costs = [-∫(x(t)-u(t)), -x(tf)]
    consolidate(u) = u[1] + u[2]
    @named rocket = ODESystem(eqs, t; costs, consolidate)
    rocket, input_idxs = structural_simplify(rocket, ([u(t)], []))

    u0map = [x(t) => 17.5]
    pmap = [u(t) => 0., tf => 8]
    jprob = JuMPControlProblem(rocket, u0map, (0, tf), pmap; steps = 201)
    jsol = solve(jprob, Ipopt.Optimizer, :Tsitouras5)
    @test isapprox(jsol.sol.t[end], 10.0, rtol = 1e-3)

    iprob = InfiniteOptControlProblem(rocket, u0map, (0, tf), pmap; steps = 200)
    isol = solve(iprob, Ipopt.Optimizer)
    @test isapprox(isol.sol.t[end], 10.0, rtol = 1e-3)
end

#@testset "Constrained optimal control problems" begin
#end
