using ModelingToolkit
import InfiniteOpt
using DiffEqDevTools, DiffEqBase
using SimpleDiffEq
using OrdinaryDiffEqSDIRK, OrdinaryDiffEqVerner, OrdinaryDiffEqTsit5, OrdinaryDiffEqFIRK
using Ipopt
using DataInterpolations
using CasADi
using Pyomo

import DiffEqBase: solve
const M = ModelingToolkit

const ENABLE_CASADI = VERSION >= v"1.11"

@testset "ODE Solution, no cost" begin
    # Test solving without anything attached.
    @parameters α=1.5 β=1.0 γ=3.0 δ=1.0
    @variables x(..) y(..)
    t = M.t_nounits
    D = M.D_nounits

    eqs = [D(x(t)) ~ α * x(t) - β * x(t) * y(t),
        D(y(t)) ~ -γ * y(t) + δ * x(t) * y(t)]

    @mtkcompile sys = System(eqs, t)
    tspan = (0.0, 1.0)
    u0map = [x(t) => 4.0, y(t) => 2.0]
    parammap = [α => 1.5, β => 1.0, γ => 3.0, δ => 1.0]

    # Test explicit method.
    jprob = JuMPDynamicOptProblem(sys, [u0map; parammap], tspan, dt = 0.01)
    jsol = solve(jprob, JuMPCollocation(Ipopt.Optimizer, constructRK4()))
    oprob = ODEProblem(sys, [u0map; parammap], tspan)
    osol = solve(oprob, SimpleRK4(), dt = 0.01)

    @test jsol.sol.u ≈ osol.u
    if ENABLE_CASADI
        cprob = CasADiDynamicOptProblem(sys, [u0map; parammap], tspan, dt = 0.01)
        csol = solve(cprob, CasADiCollocation("ipopt", constructRK4()))
        @test csol.sol.u ≈ osol.u
    end

    # Implicit method.
    osol2 = solve(oprob, ImplicitEuler(), dt = 0.01, adaptive = false)
    jsol2 = solve(jprob, JuMPCollocation(Ipopt.Optimizer, constructImplicitEuler()))
    @test ≈(jsol2.sol.u, osol2.u, rtol = 0.001)
    iprob = InfiniteOptDynamicOptProblem(sys, [u0map; parammap], tspan, dt = 0.01)
    isol = solve(iprob, InfiniteOptCollocation(Ipopt.Optimizer))
    @test ≈(isol.sol.u, osol2.u, rtol = 0.001)
    if ENABLE_CASADI
        csol2 = solve(cprob, CasADiCollocation("ipopt", constructImplicitEuler()))
        @test ≈(csol2.sol.u, osol2.u, rtol = 0.001)
    end
    pprob = PyomoDynamicOptProblem(sys, [u0map; parammap], tspan, dt = 0.01)
    psol = solve(pprob, PyomoCollocation("ipopt", BackwardEuler()))
    @test all([≈(psol.sol(t), osol2(t), rtol = 1e-2) for t in 0.0:0.01:1.0])

    # With a constraint
    u0map = Pair[]
    guess = [x(t) => 4.0, y(t) => 2.0]
    constr = [x(0.6) ~ 3.5, x(0.3) ~ 7.0]
    @mtkcompile lksys = System(eqs, t; constraints = constr)

    jprob = JuMPDynamicOptProblem(
        lksys, [u0map; parammap], tspan; guesses = guess, dt = 0.01)
    jsol = solve(jprob, JuMPCollocation(Ipopt.Optimizer, constructTsitouras5()))
    @test jsol.sol(0.6; idxs = x(t)) ≈ 3.5
    @test jsol.sol(0.3; idxs = x(t)) ≈ 7.0

    if ENABLE_CASADI
        cprob = CasADiDynamicOptProblem(
            lksys, [u0map; parammap], tspan; guesses = guess, dt = 0.01)
        csol = solve(cprob, CasADiCollocation("ipopt", constructTsitouras5()))
        @test csol.sol(0.6; idxs = x(t)) ≈ 3.5
        @test csol.sol(0.3; idxs = x(t)) ≈ 7.0
    end

    pprob = PyomoDynamicOptProblem(
        lksys, [u0map; parammap], tspan; guesses = guess, dt = 0.01)
    psol = solve(pprob, PyomoCollocation("ipopt", LagrangeLegendre(3)))
    @test psol.sol(0.6; idxs = x(t)) ≈ 3.5
    @test psol.sol(0.3; idxs = x(t)) ≈ 7.0

    iprob = InfiniteOptDynamicOptProblem(
        lksys, [u0map; parammap], tspan; guesses = guess, dt = 0.01)
    isol = solve(iprob,
        InfiniteOptCollocation(Ipopt.Optimizer, InfiniteOpt.OrthogonalCollocation(3))) # 48.564 ms, 9.58 MiB
    sol = isol.sol
    @test sol(0.6; idxs = x(t)) ≈ 3.5
    @test sol(0.3; idxs = x(t)) ≈ 7.0

    # Test whole-interval constraints
    constr = [x(t) ≳ 1, y(t) ≳ 1]
    @mtkcompile lksys = System(eqs, t; constraints = constr)
    iprob = InfiniteOptDynamicOptProblem(
        lksys, [u0map; parammap], tspan; guesses = guess, dt = 0.01)
    isol = solve(iprob,
        InfiniteOptCollocation(Ipopt.Optimizer, InfiniteOpt.OrthogonalCollocation(3)))
    @test all(u -> u > [1, 1], isol.sol.u)

    jprob = JuMPDynamicOptProblem(
        lksys, [u0map; parammap], tspan; guesses = guess, dt = 0.01)
    jsol = solve(jprob, JuMPCollocation(Ipopt.Optimizer, constructRadauIA3()))
    @test all(u -> u > [1, 1], jsol.sol.u)

    pprob = PyomoDynamicOptProblem(
        lksys, [u0map; parammap], tspan; guesses = guess, dt = 0.01)
    psol = solve(pprob, PyomoCollocation("ipopt", MidpointEuler()))
    @test all(u -> u > [1, 1], psol.sol.u)

    if ENABLE_CASADI
        cprob = CasADiDynamicOptProblem(
            lksys, [u0map; parammap], tspan; guesses = guess, dt = 0.01)
        csol = solve(cprob, CasADiCollocation("ipopt", constructRadauIA3()))
        @test all(u -> u > [1, 1], csol.sol.u)
    end
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
    @variables x(..) v(..)
    @variables u(..) [bounds = (-1.0, 1.0), input = true]
    constr = [v(1.0) ~ 0.0]
    cost = [-x(1.0)] # Maximize the final distance.
    @named block = System(
        [D(x(t)) ~ v(t), D(v(t)) ~ u(t)], t; costs = cost, constraints = constr)
    block = mtkcompile(block; inputs = [u(t)])

    u0map = [x(t) => 0.0, v(t) => 0.0]
    tspan = (0.0, 1.0)
    parammap = [u(t) => 0.0]
    jprob = JuMPDynamicOptProblem(block, [u0map; parammap], tspan; dt = 0.01)
    jsol = solve(
        jprob, JuMPCollocation(Ipopt.Optimizer, constructVerner8()), verbose = true)
    # Linear systems have bang-bang controls
    @test is_bangbang(jsol.input_sol, [-1.0], [1.0])
    # Test reached final position.
    @test ≈(jsol.sol[x(t)][end], 0.25, rtol = 1e-5)

    if ENABLE_CASADI
        cprob = CasADiDynamicOptProblem(block, [u0map; parammap], tspan; dt = 0.01)
        csol = solve(cprob, CasADiCollocation("ipopt", constructVerner8()))
        @test is_bangbang(csol.input_sol, [-1.0], [1.0])
        # Test reached final position.
        @test ≈(csol.sol[x(t)][end], 0.25, rtol = 1e-5)
    end

    # Test dynamics
    @parameters (u_interp::ConstantInterpolation)(..)
    @mtkcompile block_ode = System([D(x(t)) ~ v(t), D(v(t)) ~ u_interp(t)], t)
    spline = ctrl_to_spline(jsol.input_sol, ConstantInterpolation)
    oprob = ODEProblem(block_ode, [u0map; [u_interp => spline]], tspan)
    osol = solve(oprob, Vern8(), dt = 0.01, adaptive = false)
    @test ≈(jsol.sol.u, osol.u, rtol = 0.05)
    if ENABLE_CASADI
        @test ≈(csol.sol.u, osol.u, rtol = 0.05)
    end

    iprob = InfiniteOptDynamicOptProblem(block, [u0map; parammap], tspan; dt = 0.01)
    isol = solve(iprob, InfiniteOptCollocation(Ipopt.Optimizer))
    @test is_bangbang(isol.input_sol, [-1.0], [1.0])
    @test ≈(isol.sol[x(t)][end], 0.25, rtol = 1e-5)

    pprob = PyomoDynamicOptProblem(block, [u0map; parammap], tspan; dt = 0.01)
    psol = solve(pprob, PyomoCollocation("ipopt", BackwardEuler()))
    @test is_bangbang(psol.input_sol, [-1.0], [1.0])
    @test ≈(psol.sol[x(t)][end], 0.25, rtol = 1e-3)

    spline = ctrl_to_spline(isol.input_sol, ConstantInterpolation)
    oprob = ODEProblem(block_ode, [u0map; u_interp => spline], tspan)
    @test ≈(isol.sol.u, osol.u, rtol = 0.05)
    @test all([≈(psol.sol(t), osol(t), rtol = 0.05) for t in 0.0:0.01:1.0])

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

    @named beesys = System(eqs, t; costs)
    beesys = mtkcompile(beesys; inputs = [α])
    u0map = [w(t) => 40, q(t) => 2]
    pmap = [b => 1, c => 1, μ => 1, s => 1, ν => 1, α => 1]

    jprob = JuMPDynamicOptProblem(beesys, [u0map; pmap], tspan, dt = 0.01)
    jsol = solve(jprob, JuMPCollocation(Ipopt.Optimizer, constructTsitouras5()))
    @test is_bangbang(jsol.input_sol, [0.0], [1.0])
    iprob = InfiniteOptDynamicOptProblem(beesys, [u0map; pmap], tspan, dt = 0.01)
    isol = solve(iprob, InfiniteOptCollocation(Ipopt.Optimizer))
    @test is_bangbang(isol.input_sol, [0.0], [1.0])
    if ENABLE_CASADI
        cprob = CasADiDynamicOptProblem(beesys, [u0map; pmap], tspan; dt = 0.01)
        csol = solve(cprob, CasADiCollocation("ipopt", constructTsitouras5()))
        @test is_bangbang(csol.input_sol, [0.0], [1.0])
    end
    pprob = PyomoDynamicOptProblem(beesys, [u0map; pmap], tspan, dt = 0.01)
    psol = solve(pprob, PyomoCollocation("ipopt", BackwardEuler()))
    @test is_bangbang(psol.input_sol, [0.0], [1.0])

    @parameters (α_interp::LinearInterpolation)(..)
    eqs = [D(w(t)) ~ -μ * w(t) + b * s * α_interp(t) * w(t),
        D(q(t)) ~ -ν * q(t) + c * (1 - α_interp(t)) * s * w(t)]
    @mtkcompile beesys_ode = System(eqs, t)
    oprob = ODEProblem(beesys_ode,
        merge(Dict(u0map), Dict(pmap),
            Dict(α_interp => ctrl_to_spline(jsol.input_sol, LinearInterpolation))),
        tspan)
    osol = solve(oprob, Tsit5(); dt = 0.01, adaptive = false)
    @test ≈(osol.u, jsol.sol.u, rtol = 0.01)
    if ENABLE_CASADI
        @test ≈(osol.u, csol.sol.u, rtol = 0.01)
    end
    osol2 = solve(oprob, ImplicitEuler(); dt = 0.01, adaptive = false)
    @test ≈(osol2.u, isol.sol.u, rtol = 0.01)
    @test all([≈(psol.sol(t), osol2(t), rtol = 0.01) for t in 0.0:0.01:4.0])
end

@testset "Rocket launch" begin
    t = M.t_nounits
    D = M.D_nounits

    @parameters h_c m₀ h₀ g₀ D_c c Tₘ m_c
    @variables h(..) v(..) m(..) = m₀ [bounds = (m_c, 1)] T(..) [input = true, bounds = (0, Tₘ)]
    drag(h, v) = D_c * v^2 * exp(-h_c * (h - h₀) / h₀)
    gravity(h) = g₀ * (h₀ / h)

    eqs = [D(h(t)) ~ v(t),
        D(v(t)) ~ (T(t) - drag(h(t), v(t))) / m(t) - gravity(h(t)),
        D(m(t)) ~ -T(t) / c]

    (ts, te) = (0.0, 0.2)
    costs = [-h(te)]
    cons = [T(te) ~ 0, m(te) ~ m_c]
    @named rocket = System(eqs, t; costs, constraints = cons)
    rocket = mtkcompile(rocket; inputs = [T(t)])

    u0map = [h(t) => h₀, m(t) => m₀, v(t) => 0]
    pmap = [
        g₀ => 1, m₀ => 1.0, h_c => 500, c => 0.5 * √(g₀ * h₀), D_c => 0.5 * 620 * m₀ / g₀,
        Tₘ => 3.5 * g₀ * m₀, T(t) => 0.0, h₀ => 1, m_c => 0.6]
    jprob = JuMPDynamicOptProblem(rocket, [u0map; pmap], (ts, te); dt = 0.001, cse = false)
    jsol = solve(jprob, JuMPCollocation(Ipopt.Optimizer, constructRadauIIA5()))
    @test jsol.sol[h(t)][end] > 1.012

    if ENABLE_CASADI
        cprob = CasADiDynamicOptProblem(
            rocket, [u0map; pmap], (ts, te); dt = 0.001, cse = false)
        csol = solve(cprob, CasADiCollocation("ipopt"))
        @test csol.sol[h(t)][end] > 1.012
    end

    iprob = InfiniteOptDynamicOptProblem(rocket, [u0map; pmap], (ts, te); dt = 0.001)
    isol = solve(iprob, InfiniteOptCollocation(Ipopt.Optimizer))
    @test isol.sol[h(t)][end] > 1.012

    pprob = PyomoDynamicOptProblem(rocket, [u0map; pmap], (ts, te); dt = 0.001, cse = false)
    psol = solve(pprob, PyomoCollocation("ipopt", LagrangeRadau(4)))
    @test psol.sol[h(t)][end] > 1.012

    # Test solution
    @parameters (T_interp::CubicSpline)(..)
    eqs = [D(h(t)) ~ v(t),
        D(v(t)) ~ (T_interp(t) - drag(h(t), v(t))) / m(t) - gravity(h(t)),
        D(m(t)) ~ -T_interp(t) / c]
    @mtkcompile rocket_ode = System(eqs, t)
    interpmap = Dict(T_interp => ctrl_to_spline(jsol.input_sol, CubicSpline))
    oprob = ODEProblem(rocket_ode, merge(Dict(u0map), Dict(pmap), interpmap), (ts, te))
    osol = solve(oprob, RadauIIA5(); adaptive = false, dt = 0.001)
    @test ≈(jsol.sol.u, osol.u, rtol = 0.02)
    if ENABLE_CASADI
        @test ≈(csol.sol.u, osol.u, rtol = 0.02)
    end

    interpmap1 = Dict(T_interp => ctrl_to_spline(isol.input_sol, CubicSpline))
    oprob1 = ODEProblem(rocket_ode, merge(Dict(u0map), Dict(pmap), interpmap1), (ts, te))
    osol1 = solve(oprob1, ImplicitEuler(); adaptive = false, dt = 0.001)
    @test ≈(isol.sol.u, osol1.u, rtol = 0.01)

    interpmap2 = Dict(T_interp => ctrl_to_spline(psol.input_sol, CubicSpline))
    oprob2 = ODEProblem(rocket_ode, merge(Dict(u0map), Dict(pmap), interpmap2), (ts, te))
    osol2 = solve(oprob2, RadauIIA5(); adaptive = false, dt = 0.001)
    @test all([≈(psol.sol(t), osol2(t), rtol = 0.01) for t in 0:0.001:0.2])
end

@testset "Free final time problems" begin
    t = M.t_nounits
    D = M.D_nounits

    @variables x(..) u(..) [input = true, bounds = (0, 1)]
    @parameters tf
    eqs = [D(x(t)) ~ -2 + 0.5 * u(t)]
    # Integral cost function
    costs = [-Symbolics.Integral(t in (0, tf))(x(t) - u(t)), -x(tf)]
    consolidate(u, sub) = u[1] + u[2] + sum(sub)
    @named rocket = System(eqs, t; costs, consolidate)
    rocket = mtkcompile(rocket; inputs = [u(t)])

    u0map = [x(t) => 17.5]
    pmap = [u(t) => 0.0, tf => 8]
    jprob = JuMPDynamicOptProblem(rocket, [u0map; pmap], (0, tf); steps = 201)
    jsol = solve(jprob, JuMPCollocation(Ipopt.Optimizer, constructTsitouras5()))
    @test isapprox(jsol.sol.t[end], 10.0, rtol = 1e-3)
    @test ≈(M.objective_value(jsol), -92.75, atol = 0.25)

    if ENABLE_CASADI
        cprob = CasADiDynamicOptProblem(rocket, [u0map; pmap], (0, tf); steps = 201)
        csol = solve(cprob, CasADiCollocation("ipopt", constructTsitouras5()))
        @test isapprox(csol.sol.t[end], 10.0, rtol = 1e-3)
        @test ≈(M.objective_value(csol), -92.75, atol = 0.25)
    end

    iprob = InfiniteOptDynamicOptProblem(rocket, [u0map; pmap], (0, tf); steps = 200)
    isol = solve(iprob, InfiniteOptCollocation(Ipopt.Optimizer))
    @test isapprox(isol.sol.t[end], 10.0, rtol = 1e-3)
    @test ≈(M.objective_value(isol), -92.75, atol = 0.25)

    pprob = PyomoDynamicOptProblem(rocket, [u0map; pmap], (0, tf); steps = 201)
    psol = solve(pprob, PyomoCollocation("ipopt", BackwardEuler()))
    @test isapprox(psol.sol.t[end], 10.0, rtol = 1e-3)
    @test ≈(M.objective_value(psol), -92.75, atol = 0.1)

    @variables x(..) v(..)
    @variables u(..) [bounds = (-1.0, 1.0), input = true]
    @parameters tf
    constr = [v(tf) ~ 0, x(tf) ~ 0]
    cost = [tf] # Minimize time

    @named block = System(
        [D(x(t)) ~ v(t), D(v(t)) ~ u(t)], t; costs = cost, constraints = constr)
    block = mtkcompile(block, inputs = [u(t)])

    u0map = [x(t) => 1.0, v(t) => 0.0]
    tspan = (0.0, tf)
    parammap = [u(t) => 1.0, tf => 1.0]
    jprob = JuMPDynamicOptProblem(block, [u0map; parammap], tspan; steps = 51)
    jsol = solve(jprob, JuMPCollocation(Ipopt.Optimizer, constructVerner8()))
    @test isapprox(jsol.sol.t[end], 2.0, atol = 1e-5)

    if ENABLE_CASADI
        cprob = CasADiDynamicOptProblem(block, [u0map; parammap], (0, tf); steps = 51)
        csol = solve(cprob, CasADiCollocation("ipopt", constructVerner8()))
        @test isapprox(csol.sol.t[end], 2.0, atol = 1e-5)
    end

    iprob = InfiniteOptDynamicOptProblem(block, [u0map; parammap], tspan; steps = 51)
    isol = solve(iprob, InfiniteOptCollocation(Ipopt.Optimizer), verbose = true)
    @test isapprox(isol.sol.t[end], 2.0, atol = 1e-5)

    pprob = PyomoDynamicOptProblem(block, [u0map; parammap], tspan; steps = 51)
    psol = solve(pprob, PyomoCollocation("ipopt", BackwardEuler()))
    @test isapprox(psol.sol.t[end], 2.0, atol = 1e-5)
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

    @named cartpole = System(eqs, t; costs, constraints = cons)
    cartpole = mtkcompile(cartpole; inputs = [u])

    u0map = [D(x(t)) => 0.0, D(θ(t)) => 0.0, θ(t) => 0.0, x(t) => 0.0]
    pmap = [mₖ => 1.0, mₚ => 0.2, l => 0.5, g => 9.81, u => 0]
    jprob = JuMPDynamicOptProblem(cartpole, [u0map; pmap], tspan; dt = 0.04)
    jsol = solve(jprob, JuMPCollocation(Ipopt.Optimizer, constructRK4()))
    @test jsol.sol.u[end] ≈ [π, 0, 0, 0]

    if ENABLE_CASADI
        cprob = CasADiDynamicOptProblem(cartpole, [u0map; pmap], tspan; dt = 0.04)
        csol = solve(cprob, CasADiCollocation("ipopt", constructRK4()))
        @test csol.sol.u[end] ≈ [π, 0, 0, 0]
    end

    iprob = InfiniteOptDynamicOptProblem(cartpole, [u0map; pmap], tspan; dt = 0.04)
    isol = solve(iprob,
        InfiniteOptCollocation(Ipopt.Optimizer, InfiniteOpt.OrthogonalCollocation(2)))
    @test isol.sol.u[end] ≈ [π, 0, 0, 0]

    pprob = PyomoDynamicOptProblem(cartpole, [u0map; pmap], tspan; dt = 0.04)
    psol = solve(pprob, PyomoCollocation("ipopt", LagrangeLegendre(4)))
    @test psol.sol.u[end] ≈ [π, 0, 0, 0]
end
