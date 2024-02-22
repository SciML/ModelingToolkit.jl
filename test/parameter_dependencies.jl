using ModelingToolkit
using Test
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq
using StochasticDiffEq
using JumpProcesses
using StableRNGs
using SciMLStructures: canonicalize, Tunable, replace, replace!
using SymbolicIndexingInterface

@testset "ODESystem with callbacks" begin
    @parameters p1=1.0 p2=1.0
    @variables x(t)
    cb1 = [x ~ 2.0] => [p1 ~ 2.0] # triggers at t=-2+√6
    function affect1!(integ, u, p, ctx)
        integ.ps[p[1]] = integ.ps[p[2]]
    end
    cb2 = [x ~ 4.0] => (affect1!, [], [p1, p2], [p1]) # triggers at t=-2+√7
    cb3 = [1.0] => [p1 ~ 5.0]

    @mtkbuild sys = ODESystem(
        [D(x) ~ p1 * t + p2],
        t;
        parameter_dependencies = [p2 => 2p1],
        continuous_events = [cb1, cb2],
        discrete_events = [cb3]
    )
    @test isequal(only(parameters(sys)), p1)
    @test Set(full_parameters(sys)) == Set([p1, p2])
    prob = ODEProblem(sys, [x => 1.0], (0.0, 1.5), jac = true)
    @test prob.ps[p1] == 1.0
    @test prob.ps[p2] == 2.0
    @test_nowarn solve(prob, Tsit5())
    prob = ODEProblem(sys, [x => 1.0], (0.0, 1.5), [p1 => 1.0], jac = true)
    @test prob.ps[p1] == 1.0
    @test prob.ps[p2] == 2.0
    integ = init(prob, Tsit5())
    @test integ.ps[p1] == 1.0
    @test integ.ps[p2] == 2.0
    step!(integ, 0.5, true) # after cb1, before cb2
    @test integ.ps[p1] == 2.0
    @test integ.ps[p2] == 4.0
    step!(integ, 0.4, true) # after cb2, before cb3
    @test integ.ps[p1] == 4.0
    @test integ.ps[p2] == 8.0
    step!(integ, 0.2, true) # after cb3
    @test integ.ps[p1] == 5.0
    @test integ.ps[p2] == 10.0
end

@testset "Clock system" begin
    dt = 0.1
    @variables x(t) y(t) u(t) yd(t) ud(t) r(t) z(t)
    @parameters kp kq
    d = Clock(t, dt)
    k = ShiftIndex(d)

    eqs = [yd ~ Sample(t, dt)(y)
           ud ~ kp * (r - yd) + kq * z
           r ~ 1.0
           u ~ Hold(ud)
           D(x) ~ -x + u
           y ~ x
           z(k + 2) ~ z(k) + yd]
    @mtkbuild sys = ODESystem(eqs, t; parameter_dependencies = [kq => 2kp])

    Tf = 1.0
    prob = ODEProblem(sys, [x => 0.0, y => 0.0], (0.0, Tf),
        [kp => 1.0; z => 3.0; z(k + 1) => 2.0])
    @test_nowarn solve(prob, Tsit5(); kwargshandle = KeywordArgSilent)

    @mtkbuild sys = ODESystem(eqs, t; parameter_dependencies = [kq => 2kp],
        discrete_events = [[0.5] => [kp ~ 2.0]])
    prob = ODEProblem(sys, [x => 0.0, y => 0.0], (0.0, Tf),
        [kp => 1.0; z => 3.0; z(k + 1) => 2.0])
    @test prob.ps[kp] == 1.0
    @test prob.ps[kq] == 2.0
    @test_nowarn solve(prob, Tsit5(), kwargshandle = KeywordArgSilent)
    prob = ODEProblem(sys, [x => 0.0, y => 0.0], (0.0, Tf),
        [kp => 1.0; z => 3.0; z(k + 1) => 2.0])
    integ = init(prob, Tsit5(), kwargshandle = KeywordArgSilent)
    @test integ.ps[kp] == 1.0
    @test integ.ps[kq] == 2.0
    step!(integ, 0.6)
    @test integ.ps[kp] == 2.0
    @test integ.ps[kq] == 4.0
end

@testset "SDESystem" begin
    @parameters σ ρ β
    @variables x(t) y(t) z(t)

    eqs = [D(x) ~ σ * (y - x),
        D(y) ~ x * (ρ - z) - y,
        D(z) ~ x * y - β * z]

    noiseeqs = [0.1 * x,
        0.1 * y,
        0.1 * z]

    @named sys = ODESystem(eqs, t)
    @named sdesys = SDESystem(sys, noiseeqs; parameter_dependencies = [ρ => 2σ])
    sdesys = complete(sdesys)
    @test Set(parameters(sdesys)) == Set([σ, β])
    @test Set(full_parameters(sdesys)) == Set([σ, β, ρ])

    prob = SDEProblem(
        sdesys, [x => 1.0, y => 0.0, z => 0.0], (0.0, 100.0), [σ => 10.0, β => 2.33])
    @test prob.ps[ρ] == 2prob.ps[σ]
    @test_nowarn solve(prob, SRIW1())

    @named sys = ODESystem(eqs, t)
    @named sdesys = SDESystem(sys, noiseeqs; parameter_dependencies = [ρ => 2σ],
        discrete_events = [[10.0] => [σ ~ 15.0]])
    sdesys = complete(sdesys)
    prob = SDEProblem(
        sdesys, [x => 1.0, y => 0.0, z => 0.0], (0.0, 100.0), [σ => 10.0, β => 2.33])
    integ = init(prob, SRIW1())
    @test integ.ps[σ] == 10.0
    @test integ.ps[ρ] == 20.0
    step!(integ, 11.0)
    @test integ.ps[σ] == 15.0
    @test integ.ps[ρ] == 30.0
end

@testset "JumpSystem" begin
    rng = StableRNG(12345)
    @parameters β γ
    @constants h = 1
    @variables S(t) I(t) R(t)
    rate₁ = β * S * I * h
    affect₁ = [S ~ S - 1 * h, I ~ I + 1]
    rate₃ = γ * I * h
    affect₃ = [I ~ I * h - 1, R ~ R + 1]
    j₁ = ConstantRateJump(rate₁, affect₁)
    j₃ = ConstantRateJump(rate₃, affect₃)
    @named js2 = JumpSystem(
        [j₁, j₃], t, [S, I, R], [γ]; parameter_dependencies = [β => 0.01γ])
    @test isequal(only(parameters(js2)), γ)
    @test Set(full_parameters(js2)) == Set([γ, β])
    js2 = complete(js2)
    tspan = (0.0, 250.0)
    u₀map = [S => 999, I => 1, R => 0]
    parammap = [γ => 0.01]
    dprob = DiscreteProblem(js2, u₀map, tspan, parammap)
    jprob = JumpProblem(js2, dprob, Direct(), save_positions = (false, false), rng = rng)
    @test jprob.ps[γ] == 0.01
    @test jprob.ps[β] == 0.0001
    @test_nowarn solve(jprob, SSAStepper())

    @named js2 = JumpSystem(
        [j₁, j₃], t, [S, I, R], [γ]; parameter_dependencies = [β => 0.01γ],
        discrete_events = [[10.0] => [γ ~ 0.02]])
    js2 = complete(js2)
    dprob = DiscreteProblem(js2, u₀map, tspan, parammap)
    jprob = JumpProblem(js2, dprob, Direct(), save_positions = (false, false), rng = rng)
    integ = init(jprob, SSAStepper())
    @test integ.ps[γ] == 0.01
    @test integ.ps[β] == 0.0001
    step!(integ, 11.0)
    @test integ.ps[γ] == 0.02
    @test integ.ps[β] == 0.0002
end

@testset "SciMLStructures interface" begin
    @parameters p1=1.0 p2=1.0
    @variables x(t)
    cb1 = [x ~ 2.0] => [p1 ~ 2.0] # triggers at t=-2+√6
    function affect1!(integ, u, p, ctx)
        integ.ps[p[1]] = integ.ps[p[2]]
    end
    cb2 = [x ~ 4.0] => (affect1!, [], [p1, p2], [p1]) # triggers at t=-2+√7
    cb3 = [1.0] => [p1 ~ 5.0]

    @mtkbuild sys = ODESystem(
        [D(x) ~ p1 * t + p2],
        t;
        parameter_dependencies = [p2 => 2p1]
    )
    prob = ODEProblem(sys, [x => 1.0], (0.0, 1.5), [p1 => 1.0], jac = true)
    prob.ps[p1] = 3.0
    @test prob.ps[p1] == 3.0
    @test prob.ps[p2] == 6.0

    ps = prob.p
    buffer, repack, _ = canonicalize(Tunable(), ps)
    @test only(buffer) == 3.0
    buffer[1] = 4.0
    ps = repack(buffer)
    @test getp(sys, p1)(ps) == 4.0
    @test getp(sys, p2)(ps) == 8.0

    replace!(Tunable(), ps, [1.0])
    @test getp(sys, p1)(ps) == 1.0
    @test getp(sys, p2)(ps) == 2.0

    ps2 = replace(Tunable(), ps, [2.0])
    @test getp(sys, p1)(ps2) == 2.0
    @test getp(sys, p2)(ps2) == 4.0
end
