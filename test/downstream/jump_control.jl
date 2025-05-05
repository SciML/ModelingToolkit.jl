using ModelingToolkit
import JuMP, InfiniteOpt
using DiffEqDevTools, DiffEqBase
using SimpleDiffEq
using OrdinaryDiffEqSDIRK, OrdinaryDiffEqVerner, OrdinaryDiffEqTsit5, OrdinaryDiffEqFIRK
using Ipopt
using DataInterpolations
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
    jprob = JuMPDynamicOptProblem(sys, u0map, tspan, parammap, dt = 0.01)
    @test JuMP.num_constraints(jprob.model) == 2 # initials
    jsol = solve(jprob, Ipopt.Optimizer, constructRK4, silent = true)
    oprob = ODEProblem(sys, u0map, tspan, parammap)
    osol = solve(oprob, SimpleRK4(), dt = 0.01)
    @test jsol.sol.u ≈ osol.u

    # Implicit method.
    jsol2 = solve(jprob, Ipopt.Optimizer, constructImplicitEuler, silent = true) # 63.031 ms, 26.49 MiB
    osol2 = solve(oprob, ImplicitEuler(), dt = 0.01, adaptive = false) # 129.375 μs, 61.91 KiB
    @test ≈(jsol2.sol.u, osol2.u, rtol = 0.001)
    iprob = InfiniteOptDynamicOptProblem(sys, u0map, tspan, parammap, dt = 0.01)
    isol = solve(iprob, Ipopt.Optimizer,
        derivative_method = InfiniteOpt.FiniteDifference(InfiniteOpt.Backward()),
        silent = true) # 11.540 ms, 4.00 MiB
    @test ≈(jsol2.sol.u, osol2.u, rtol = 0.001)

    # With a constraint
    u0map = Pair[]
    guess = [x(t) => 4.0, y(t) => 2.0]
    constr = [x(0.6) ~ 3.5, x(0.3) ~ 7.0]
    @mtkbuild lksys = ODESystem(eqs, t; constraints = constr)

    jprob = JuMPDynamicOptProblem(lksys, u0map, tspan, parammap; guesses = guess, dt = 0.01)
    @test JuMP.num_constraints(jprob.model) == 2
    jsol = solve(jprob, Ipopt.Optimizer, constructTsitouras5, silent = true) # 12.190 s, 9.68 GiB
    @test jsol.sol(0.6)[1] ≈ 3.5
    @test jsol.sol(0.3)[1] ≈ 7.0

    iprob = InfiniteOptDynamicOptProblem(
        lksys, u0map, tspan, parammap; guesses = guess, dt = 0.01)
    isol = solve(iprob, Ipopt.Optimizer,
        derivative_method = InfiniteOpt.OrthogonalCollocation(3), silent = true) # 48.564 ms, 9.58 MiB
    sol = isol.sol
    @test sol(0.6)[1] ≈ 3.5
    @test sol(0.3)[1] ≈ 7.0

    # Test whole-interval constraints
    constr = [x(t) ≳ 1, y(t) ≳ 1]
    @mtkbuild lksys = ODESystem(eqs, t; constraints = constr)
    iprob = InfiniteOptDynamicOptProblem(
        lksys, u0map, tspan, parammap; guesses = guess, dt = 0.01)
    isol = solve(iprob, Ipopt.Optimizer,
        derivative_method = InfiniteOpt.OrthogonalCollocation(3), silent = true) # 48.564 ms, 9.58 MiB
    @test all(u -> u > [1, 1], isol.sol.u)

    jprob = JuMPDynamicOptProblem(lksys, u0map, tspan, parammap; guesses = guess, dt = 0.01)
    jsol = solve(jprob, Ipopt.Optimizer, constructRadauIA3, silent = true) # 12.190 s, 9.68 GiB
    @test all(u -> u > [1, 1], jsol.sol.u)
end

function is_bangbang(input_sol, lbounds, ubounds, rtol = 1e-4)
    for v in 1:(length(input_sol.u[1]) - 1)
        all(i -> ≈(i[v], bounds[v]; rtol) || ≈(i[v], bounds[u]; rtol), input_sol.u) ||
            return false
    end
    true
end

function ctrl_to_spline(inputsol, splineType)
    us = reduce(vcat, inputsol.u)
    ts = reduce(vcat, inputsol.t)
    splineType(us, ts)
end

@testset "Linear systems" begin
    # Double integrator
    t = M.t_nounits
    D = M.D_nounits
    @variables x(..) [bounds = (0.0, 0.25)] v(..)
    @variables u(..) [bounds = (-1.0, 1.0), input = true]
    constr = [v(1.0) ~ 0.0]
    cost = [-x(1.0)] # Maximize the final distance.
    @named block = ODESystem(
        [D(x(t)) ~ v(t), D(v(t)) ~ u(t)], t; costs = cost, constraints = constr)
    block, input_idxs = structural_simplify(block, ([u(t)], []))

    u0map = [x(t) => 0.0, v(t) => 0.0]
    tspan = (0.0, 1.0)
    parammap = [u(t) => 0.0]
    jprob = JuMPDynamicOptProblem(block, u0map, tspan, parammap; dt = 0.01)
    jsol = solve(jprob, Ipopt.Optimizer, constructVerner8, silent = true)
    # Linear systems have bang-bang controls
    @test is_bangbang(jsol.input_sol, [-1.0], [1.0])
    # Test reached final position.
    @test ≈(jsol.sol.u[end][2], 0.25, rtol = 1e-5)
    # Test dynamics
    @parameters (u_interp::ConstantInterpolation)(..)
    @mtkbuild block_ode = ODESystem([D(x(t)) ~ v(t), D(v(t)) ~ u_interp(t)], t)
    spline = ctrl_to_spline(jsol.input_sol, ConstantInterpolation)
    oprob = ODEProblem(block_ode, u0map, tspan, [u_interp => spline])
    osol = solve(oprob, Vern8(), dt = 0.01, adaptive = false)
    @test ≈(jsol.sol.u, osol.u, rtol = 0.05)

    iprob = InfiniteOptDynamicOptProblem(block, u0map, tspan, parammap; dt = 0.01)
    isol = solve(iprob, Ipopt.Optimizer; silent = true)
    @test is_bangbang(isol.input_sol, [-1.0], [1.0])
    @test ≈(isol.sol.u[end][2], 0.25, rtol = 1e-5)
    osol = solve(oprob, ImplicitEuler(); dt = 0.01, adaptive = false)
    @test ≈(isol.sol.u, osol.u, rtol = 0.05)

    ###################
    ### Bee example ###
    ###################
    # From Lawrence Evans' notes
    @variables w(..) q(..) α(t) [input = true, bounds = (0, 1)]
    @parameters b c μ s ν

    tspan = (0, 4)
    eqs = [D(w(t)) ~ -μ * w(t) + b * s * α * w(t),
        D(q(t)) ~ -ν * q(t) + c * (1 - α) * s * w(t)]
    costs = [-q(tspan[2])]

    @named beesys = ODESystem(eqs, t; costs)
    beesys, input_idxs = structural_simplify(beesys, ([α], []))
    u0map = [w(t) => 40, q(t) => 2]
    pmap = [b => 1, c => 1, μ => 1, s => 1, ν => 1, α => 1]

    jprob = JuMPDynamicOptProblem(beesys, u0map, tspan, pmap, dt = 0.01)
    jsol = solve(jprob, Ipopt.Optimizer, constructTsitouras5, silent = true)
    @test is_bangbang(jsol.input_sol, [0.0], [1.0])
    iprob = InfiniteOptDynamicOptProblem(beesys, u0map, tspan, pmap, dt = 0.01)
    isol = solve(iprob, Ipopt.Optimizer; silent = true)
    @test is_bangbang(isol.input_sol, [0.0], [1.0])

    @parameters (α_interp::LinearInterpolation)(..)
    eqs = [D(w(t)) ~ -μ * w(t) + b * s * α_interp(t) * w(t),
        D(q(t)) ~ -ν * q(t) + c * (1 - α_interp(t)) * s * w(t)]
    @mtkbuild beesys_ode = ODESystem(eqs, t)
    oprob = ODEProblem(beesys_ode,
        u0map,
        tspan,
        merge(Dict(pmap),
            Dict(α_interp => ctrl_to_spline(jsol.input_sol, LinearInterpolation))))
    osol = solve(oprob, Tsit5(); dt = 0.01, adaptive = false)
    @test ≈(osol.u, jsol.sol.u, rtol = 0.01)
    osol2 = solve(oprob, ImplicitEuler(); dt = 0.01, adaptive = false)
    @test ≈(osol2.u, isol.sol.u, rtol = 0.01)
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

    (ts, te) = (0.0, 0.2)
    costs = [-h(te)]
    cons = [T(te) ~ 0, m(te) ~ m_c]
    @named rocket = ODESystem(eqs, t; costs, constraints = cons)
    rocket, input_idxs = structural_simplify(rocket, ([T(t)], []))

    u0map = [h(t) => h₀, m(t) => m₀, v(t) => 0]
    pmap = [
        g₀ => 1, m₀ => 1.0, h_c => 500, c => 0.5 * √(g₀ * h₀), D_c => 0.5 * 620 * m₀ / g₀,
        Tₘ => 3.5 * g₀ * m₀, T(t) => 0.0, h₀ => 1, m_c => 0.6]
    jprob = JuMPDynamicOptProblem(rocket, u0map, (ts, te), pmap; dt = 0.001, cse = false)
    jsol = solve(jprob, Ipopt.Optimizer, constructRadauIIA5, silent = true)
    @test jsol.sol.u[end][1] > 1.012

    iprob = InfiniteOptDynamicOptProblem(
        rocket, u0map, (ts, te), pmap; dt = 0.001)
    isol = solve(iprob, Ipopt.Optimizer, silent = true)
    @test isol.sol.u[end][1] > 1.012

    # Test solution
    @parameters (T_interp::CubicSpline)(..)
    eqs = [D(h(t)) ~ v(t),
        D(v(t)) ~ (T_interp(t) - drag(h(t), v(t))) / m(t) - gravity(h(t)),
        D(m(t)) ~ -T_interp(t) / c]
    @mtkbuild rocket_ode = ODESystem(eqs, t)
    interpmap = Dict(T_interp => ctrl_to_spline(jsol.input_sol, CubicSpline))
    oprob = ODEProblem(rocket_ode, u0map, (ts, te), merge(Dict(pmap), interpmap))
    osol = solve(oprob, RadauIIA5(); adaptive = false, dt = 0.001)
    @test ≈(jsol.sol.u, osol.u, rtol = 0.02)

    interpmap1 = Dict(T_interp => ctrl_to_spline(isol.input_sol, CubicSpline))
    oprob1 = ODEProblem(rocket_ode, u0map, (ts, te), merge(Dict(pmap), interpmap1))
    osol1 = solve(oprob1, ImplicitEuler(); adaptive = false, dt = 0.001)
    @test ≈(isol.sol.u, osol1.u, rtol = 0.01)
end

@testset "Free final time problems" begin
    t = M.t_nounits
    D = M.D_nounits

    @variables x(..) u(..) [input = true, bounds = (0, 1)]
    @parameters tf
    eqs = [D(x(t)) ~ -2 + 0.5 * u(t)]
    # Integral cost function
    costs = [-Symbolics.Integral(t in (0, tf))(x(t) - u(t)), -x(tf)]
    consolidate(u) = u[1] + u[2]
    @named rocket = ODESystem(eqs, t; costs, consolidate)
    rocket, input_idxs = structural_simplify(rocket, ([u(t)], []))

    u0map = [x(t) => 17.5]
    pmap = [u(t) => 0.0, tf => 8]
    jprob = JuMPDynamicOptProblem(rocket, u0map, (0, tf), pmap; steps = 201)
    jsol = solve(jprob, Ipopt.Optimizer, constructTsitouras5, silent = true)
    @test isapprox(jsol.sol.t[end], 10.0, rtol = 1e-3)

    iprob = InfiniteOptDynamicOptProblem(rocket, u0map, (0, tf), pmap; steps = 200)
    isol = solve(iprob, Ipopt.Optimizer, silent = true)
    @test isapprox(isol.sol.t[end], 10.0, rtol = 1e-3)

    @variables x(..) v(..)
    @variables u(..) [bounds = (-1.0, 1.0), input = true]
    @parameters tf
    constr = [v(tf) ~ 0, x(tf) ~ 0]
    cost = [tf] # Minimize time

    @named block = ODESystem(
        [D(x(t)) ~ v(t), D(v(t)) ~ u(t)], t; costs = cost, constraints = constr)
    block, input_idxs = structural_simplify(block, ([u(t)], []))

    u0map = [x(t) => 1.0, v(t) => 0.0]
    tspan = (0.0, tf)
    parammap = [u(t) => 0.0, tf => 1.0]
    jprob = JuMPDynamicOptProblem(block, u0map, tspan, parammap; steps = 51)
    jsol = solve(jprob, Ipopt.Optimizer, constructVerner8, silent = true)
    @test isapprox(jsol.sol.t[end], 2.0, atol = 1e-5)

    iprob = InfiniteOptDynamicOptProblem(block, u0map, tspan, parammap; steps = 51)
    isol = solve(iprob, Ipopt.Optimizer, silent = true)
    @test isapprox(isol.sol.t[end], 2.0, atol = 1e-5)
end

@testset "Cart-pole problem" begin
    t = M.t_nounits
    D = M.D_nounits
    # gravity, length, moment of Inertia, drag coeff
    @parameters g l mₚ mₖ
    @variables x(..) θ(..) u(t) [input = true, bounds = (-10, 10)]

    s = sin(θ(t))
    c = cos(θ(t))
    H = [mₖ+mₚ mₚ*l*c
         mₚ*l*c mₚ*l^2]
    C = [0 -mₚ*D(θ(t))*l*s
         0 0]
    qd = [D(x(t)), D(θ(t))]
    G = [0, mₚ * g * l * s]
    B = [1, 0]

    tf = 5
    rhss = -H \ Vector(C * qd + G - B * u)
    eqs = [D(D(x(t))) ~ rhss[1], D(D(θ(t))) ~ rhss[2]]
    cons = [θ(tf) ~ π, x(tf) ~ 0, D(θ(tf)) ~ 0, D(x(tf)) ~ 0]
    costs = [Symbolics.Integral(t in (0, tf))(u^2)]
    tspan = (0, tf)

    @named cartpole = ODESystem(eqs, t; costs, constraints = cons)
    cartpole, input_idxs = structural_simplify(cartpole, ([u], []))

    u0map = [D(x(t)) => 0.0, D(θ(t)) => 0.0, θ(t) => 0.0, x(t) => 0.0]
    pmap = [mₖ => 1.0, mₚ => 0.2, l => 0.5, g => 9.81, u => 0]
    jprob = JuMPDynamicOptProblem(cartpole, u0map, tspan, pmap; dt = 0.04)
    jsol = solve(jprob, Ipopt.Optimizer, constructRK4, silent = true)
    @test jsol.sol.u[end] ≈ [π, 0, 0, 0]

    iprob = InfiniteOptDynamicOptProblem(cartpole, u0map, tspan, pmap; dt = 0.04)
    isol = solve(iprob, Ipopt.Optimizer, silent = true)
    @test isol.sol.u[end] ≈ [π, 0, 0, 0]
end
