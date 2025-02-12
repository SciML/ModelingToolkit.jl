using ModelingToolkit, DiffEqBase, JumpProcesses, Test, LinearAlgebra
using Random, StableRNGs, NonlinearSolve
using OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
MT = ModelingToolkit

rng = StableRNG(12345)

# basic MT SIR model with tweaks
@parameters β γ
@constants h = 1
@variables S(t) I(t) R(t)
rate₁ = β * S * I * h
affect₁ = [S ~ S - 1 * h, I ~ I + 1]
rate₂ = γ * I + t
affect₂ = [I ~ I - 1, R ~ R + 1]
j₁ = ConstantRateJump(rate₁, affect₁)
j₂ = VariableRateJump(rate₂, affect₂)
@named js = JumpSystem([j₁, j₂], t, [S, I, R], [β, γ])
unknowntoid = Dict(MT.value(unknown) => i for (i, unknown) in enumerate(unknowns(js)))
mtjump1 = MT.assemble_crj(js, j₁, unknowntoid)
mtjump2 = MT.assemble_vrj(js, j₂, unknowntoid)

# doc version
rate1(u, p, t) = (0.1 / 1000.0) * u[1] * u[2]
function affect1!(integrator)
    integrator.u[1] -= 1
    integrator.u[2] += 1
end
jump1 = ConstantRateJump(rate1, affect1!)
rate2(u, p, t) = 0.01u[2] + t
function affect2!(integrator)
    integrator.u[2] -= 1
    integrator.u[3] += 1
end
jump2 = VariableRateJump(rate2, affect2!)

# test crjs
u = [100, 9, 5]
p = (0.1 / 1000, 0.01)
tf = 1.0
mutable struct TestInt{U, V, T}
    u::U
    p::V
    t::T
end
mtintegrator = TestInt(u, p, tf)
integrator = TestInt(u, p, tf)
@test abs(mtjump1.rate(u, p, tf) - jump1.rate(u, p, tf)) < 10 * eps()
@test abs(mtjump2.rate(u, p, tf) - jump2.rate(u, p, tf)) < 10 * eps()
mtjump1.affect!(mtintegrator)
jump1.affect!(integrator)
@test all(integrator.u .== mtintegrator.u)
mtintegrator.u .= u;
integrator.u .= u;
mtjump2.affect!(mtintegrator)
jump2.affect!(integrator)
@test all(integrator.u .== mtintegrator.u)

# test MT can make and solve a jump problem
rate₃ = γ * I * h
affect₃ = [I ~ I * h - 1, R ~ R + 1]
j₃ = ConstantRateJump(rate₃, affect₃)
@named js2 = JumpSystem([j₁, j₃], t, [S, I, R], [β, γ])
js2 = complete(js2)
u₀ = [999, 1, 0];
p = (0.1 / 1000, 0.01);
tspan = (0.0, 250.0);
u₀map = [S => 999, I => 1, R => 0]
parammap = [β => 0.1 / 1000, γ => 0.01]
dprob = DiscreteProblem(js2, u₀map, tspan, parammap)
jprob = JumpProblem(js2, dprob, Direct(); save_positions = (false, false), rng)
Nsims = 30000
function getmean(jprob, Nsims; use_stepper = true)
    m = 0.0
    for i in 1:Nsims
        sol = use_stepper ? solve(jprob, SSAStepper()) : solve(jprob)
        m += sol[end, end]
    end
    m / Nsims
end
m = getmean(jprob, Nsims)

# test auto-alg selection works
jprobb = JumpProblem(js2, dprob; save_positions = (false, false), rng)
mb = getmean(jprobb, Nsims; use_stepper = false)
@test abs(m - mb) / m < 0.01

@variables S2(t)
obs = [S2 ~ 2 * S]
@named js2b = JumpSystem([j₁, j₃], t, [S, I, R], [β, γ], observed = obs)
js2b = complete(js2b)
dprob = DiscreteProblem(js2b, u₀map, tspan, parammap)
jprob = JumpProblem(js2b, dprob, Direct(); save_positions = (false, false), rng)
sol = solve(jprob, SSAStepper(); saveat = tspan[2] / 10)
@test all(2 .* sol[S] .== sol[S2])

# test save_positions is working
jprob = JumpProblem(js2, dprob, Direct(); save_positions = (false, false), rng)
sol = solve(jprob, SSAStepper(); saveat = 1.0)
@test all((sol.t) .== collect(0.0:tspan[2]))

#test the MT JumpProblem rates/affects are correct
rate2(u, p, t) = 0.01u[2]
jump2 = ConstantRateJump(rate2, affect2!)
mtjumps = jprob.discrete_jump_aggregation
@test abs(mtjumps.rates[1](u, (p,), tf) - jump1.rate(u, p, tf)) < 10 * eps()
@test abs(mtjumps.rates[2](u, (p,), tf) - jump2.rate(u, p, tf)) < 10 * eps()
mtjumps.affects![1](mtintegrator)
jump1.affect!(integrator)
@test all(integrator.u .== mtintegrator.u)
mtintegrator.u .= u;
integrator.u .= u;
mtjumps.affects![2](mtintegrator)
jump2.affect!(integrator)
@test all(integrator.u .== mtintegrator.u)

# direct vers
p = (0.1 / 1000, 0.01)
prob = DiscreteProblem([999, 1, 0], (0.0, 250.0), p)
r1(u, p, t) = (0.1 / 1000.0) * u[1] * u[2]
function a1!(integrator)
    integrator.u[1] -= 1
    integrator.u[2] += 1
end
j1 = ConstantRateJump(r1, a1!)
r2(u, p, t) = 0.01u[2]
function a2!(integrator)
    integrator.u[2] -= 1
    integrator.u[3] += 1
end
j2 = ConstantRateJump(r2, a2!)
jset = JumpSet((), (j1, j2), nothing, nothing)
jprob = JumpProblem(prob, Direct(), jset; save_positions = (false, false), rng)
m2 = getmean(jprob, Nsims)

# test JumpSystem solution agrees with direct version
@test abs(m - m2) / m < 0.01

# mass action jump tests for SIR model
maj1 = MassActionJump(2 * β / 2, [S => 1, I => 1], [S => -1, I => 1])
maj2 = MassActionJump(γ, [I => 1], [I => -1, R => 1])
@named js3 = JumpSystem([maj1, maj2], t, [S, I, R], [β, γ])
js3 = complete(js3)
dprob = DiscreteProblem(js3, u₀map, tspan, parammap)
jprob = JumpProblem(js3, dprob, Direct(); rng)
m3 = getmean(jprob, Nsims)
@test abs(m - m3) / m < 0.01

# maj jump test with various dep graphs
@named js3b = JumpSystem([maj1, maj2], t, [S, I, R], [β, γ])
js3b = complete(js3b)
jprobb = JumpProblem(js3b, dprob, NRM(); rng)
m4 = getmean(jprobb, Nsims)
@test abs(m - m4) / m < 0.01
jprobc = JumpProblem(js3b, dprob, RSSA(); rng)
m4 = getmean(jprobc, Nsims)
@test abs(m - m4) / m < 0.01

# mass action jump tests for other reaction types (zero order, decay)
maj1 = MassActionJump(2.0, [0 => 1], [S => 1])
maj2 = MassActionJump(γ, [S => 1], [S => -1])
@named js4 = JumpSystem([maj1, maj2], t, [S], [β, γ])
js4 = complete(js4)
dprob = DiscreteProblem(js4, [S => 999], (0, 1000.0), [β => 100.0, γ => 0.01])
jprob = JumpProblem(js4, dprob, Direct(); rng)
m4 = getmean(jprob, Nsims)
@test abs(m4 - 2.0 / 0.01) * 0.01 / 2.0 < 0.01

# test second order rx runs
maj1 = MassActionJump(2.0, [0 => 1], [S => 1])
maj2 = MassActionJump(γ, [S => 2], [S => -1])
@named js4 = JumpSystem([maj1, maj2], t, [S], [β, γ])
js4 = complete(js4)
dprob = DiscreteProblem(js4, [S => 999], (0, 1000.0), [β => 100.0, γ => 0.01])
jprob = JumpProblem(js4, dprob, Direct(); rng)
sol = solve(jprob, SSAStepper());

# issue #819
@testset "Combined system name collisions" begin
    sys1 = JumpSystem([maj1, maj2], t, [S], [β, γ], name = :sys1)
    sys2 = JumpSystem([maj1, maj2], t, [S], [β, γ], name = :sys1)
    @test_throws ArgumentError JumpSystem([sys1.γ ~ sys2.γ], t, [], [],
        systems = [sys1, sys2], name = :foo)
end

# test if param mapper is setup correctly for callbacks
let
    @parameters k1 k2 k3
    @variables A(t) B(t)
    maj1 = MassActionJump(k1 * k3, [0 => 1], [A => -1, B => 1])
    maj2 = MassActionJump(k2, [B => 1], [A => 1, B => -1])
    @named js5 = JumpSystem([maj1, maj2], t, [A, B], [k1, k2, k3])
    js5 = complete(js5)
    p = [k1 => 2.0, k2 => 0.0, k3 => 0.5]
    u₀ = [A => 100, B => 0]
    tspan = (0.0, 2000.0)
    dprob = DiscreteProblem(js5, u₀, tspan, p)
    jprob = JumpProblem(js5, dprob, Direct(); save_positions = (false, false), rng)
    @test all(jprob.massaction_jump.scaled_rates .== [1.0, 0.0])

    pcondit(u, t, integrator) = t == 1000.0
    function paffect!(integrator)
        integrator.ps[k1] = 0.0
        integrator.ps[k2] = 1.0
        reset_aggregated_jumps!(integrator)
    end
    cb = DiscreteCallback(pcondit, paffect!)
    sol = solve(jprob, SSAStepper(); tstops = [1000.0], callback = cb)
    @test sol.u[end][1] == 100
end

# observed variable handling
@variables OBS(t)
@named js5 = JumpSystem([maj1, maj2], t, [S], [β, γ]; observed = [OBS ~ 2 * S * h])
OBS2 = OBS
@test isequal(OBS2, @nonamespace js5.OBS)
@unpack OBS = js5
@test isequal(OBS2, OBS)

# test to make sure dep graphs are correct
let
    # A + 2X --> 3X
    # 3X --> A + 2X
    # B --> X
    # X --> B
    @variables A(t) X(t) B(t)
    jumps = [MassActionJump(1.0, [A => 1, X => 2], [A => -1, X => 1]),
        MassActionJump(1.0, [X => 3], [A => 1, X => -1]),
        MassActionJump(1.0, [B => 1], [B => -1, X => 1]),
        MassActionJump(1.0, [X => 1], [B => 1, X => -1])]
    @named js = JumpSystem(jumps, t, [A, X, B], [])
    jdeps = asgraph(js)
    vdeps = variable_dependencies(js)
    vtoj = jdeps.badjlist
    @test vtoj == [[1], [1, 2, 4], [3]]
    jtov = vdeps.badjlist
    @test jtov == [[1, 2], [1, 2], [2, 3], [2, 3]]
    jtoj = eqeq_dependencies(jdeps, vdeps).fadjlist
    @test jtoj == [[1, 2, 4], [1, 2, 4], [1, 2, 3, 4], [1, 2, 3, 4]]
end

# Create JumpProblems for systems without parameters
# Issue#2559
@parameters k
@variables X(t)
rate = k
affect = [X ~ X - 1]

crj = ConstantRateJump(1.0, [X ~ X - 1])
js1 = complete(JumpSystem([crj], t, [X], [k]; name = :js1))
js2 = complete(JumpSystem([crj], t, [X], []; name = :js2))

maj = MassActionJump(1.0, [X => 1], [X => -1])
js3 = complete(JumpSystem([maj], t, [X], [k]; name = :js2))
js4 = complete(JumpSystem([maj], t, [X], []; name = :js3))

u0 = [X => 10]
tspan = (0.0, 1.0)
ps = [k => 1.0]

dp1 = DiscreteProblem(js1, u0, tspan, ps)
dp2 = DiscreteProblem(js2, u0, tspan)
dp3 = DiscreteProblem(js3, u0, tspan, ps)
dp4 = DiscreteProblem(js4, u0, tspan)

@test_nowarn jp1 = JumpProblem(js1, dp1, Direct())
@test_nowarn jp2 = JumpProblem(js2, dp2, Direct())
@test_nowarn jp3 = JumpProblem(js3, dp3, Direct())
@test_nowarn jp4 = JumpProblem(js4, dp4, Direct())

# Ensure `structural_simplify` (and `@mtkbuild`) works on JumpSystem (by doing nothing)
# Issue#2558
@parameters k
@variables X(t)
rate = k
affect = [X ~ X - 1]

j1 = ConstantRateJump(k, [X ~ X - 1])
@test_nowarn @mtkbuild js1 = JumpSystem([j1], t, [X], [k])

# test correct autosolver is selected, which implies appropriate dep graphs are available
let
    @parameters k
    @variables X(t)
    rate = k
    affect = [X ~ X - 1]
    j1 = ConstantRateJump(k, [X ~ X - 1])

    Nv = [1, JumpProcesses.USE_DIRECT_THRESHOLD + 1, JumpProcesses.USE_RSSA_THRESHOLD + 1]
    algtypes = [Direct, RSSA, RSSACR]
    for (N, algtype) in zip(Nv, algtypes)
        @named jsys = JumpSystem([deepcopy(j1) for _ in 1:N], t, [X], [k])
        jsys = complete(jsys)
        dprob = DiscreteProblem(jsys, [X => 10], (0.0, 10.0), [k => 1])
        jprob = JumpProblem(jsys, dprob)
        @test jprob.aggregator isa algtype
    end
end

# basic VariableRateJump test
let
    N = 1000  # number of simulations for testing solve accuracy
    Random.seed!(rng, 1111)
    @variables A(t) B(t) C(t)
    @parameters k
    vrj = VariableRateJump(k * (sin(t) + 1), [A ~ A + 1, C ~ C + 2])
    js = complete(JumpSystem([vrj], t, [A, C], [k]; name = :js, observed = [B ~ C * A]))
    oprob = ODEProblem(js, [A => 0, C => 0], (0.0, 10.0), [k => 1.0])
    jprob = JumpProblem(js, oprob, Direct(); rng)
    sol = solve(jprob, Tsit5())

    # test observed and symbolic indexing work
    @test all(sol[:A] .* sol[:C] .== sol[:B])

    dt = 1.0
    tv = range(0.0, 10.0; step = 1.0)
    cmean = zeros(11)
    for n in 1:N
        sol = solve(jprob, Tsit5(); save_everystep = false, saveat = dt)
        cmean += Array(sol(tv; idxs = :C))
    end
    cmean ./= N

    vrjrate(u, p, t) = p[1] * (sin(t) + 1)
    function vrjaffect!(integ)
        integ.u[1] += 1
        integ.u[2] += 2
        nothing
    end
    vrj2 = VariableRateJump(vrjrate, vrjaffect!)
    oprob2 = ODEProblem((du, u, p, t) -> (du .= 0; nothing), [0, 0], (0.0, 10.0), (1.0,))
    jprob2 = JumpProblem(oprob2, Direct(), vrj2; rng)
    cmean2 = zeros(11)
    for n in 1:N
        sol2 = solve(jprob2, Tsit5(); saveat = dt)
        cmean2 += Array(sol2(tv; idxs = 2))
    end
    cmean2 ./= N

    @test all(abs.(cmean .- cmean2) .<= 0.05 .* cmean)
end

# collect_vars! tests for jumps
let
    @variables x1(t) x2(t) x3(t) x4(t) x5(t)
    @parameters p1 p2 p3 p4 p5
    j1 = ConstantRateJump(p1, [x1 ~ x1 + 1])
    j2 = MassActionJump(p2, [x2 => 1], [x3 => -1])
    j3 = VariableRateJump(p3, [x3 ~ x3 + 1, x4 ~ x4 + 1])
    j4 = MassActionJump(p4 * p5, [x1 => 1, x5 => 1], [x1 => -1, x5 => -1, x2 => 1])
    us = Set()
    ps = Set()
    iv = t

    MT.collect_vars!(us, ps, j1, iv)
    @test issetequal(us, [x1])
    @test issetequal(ps, [p1])

    empty!(us)
    empty!(ps)
    MT.collect_vars!(us, ps, j2, iv)
    @test issetequal(us, [x2, x3])
    @test issetequal(ps, [p2])

    empty!(us)
    empty!(ps)
    MT.collect_vars!(us, ps, j3, iv)
    @test issetequal(us, [x3, x4])
    @test issetequal(ps, [p3])

    empty!(us)
    empty!(ps)
    MT.collect_vars!(us, ps, j4, iv)
    @test issetequal(us, [x1, x5, x2])
    @test issetequal(ps, [p4, p5])
end

# scoping tests
let
    @variables x1(t) x2(t) x3(t) x4(t) x5(t)
    x2 = ParentScope(x2)
    x3 = ParentScope(ParentScope(x3))
    x4 = DelayParentScope(x4, 2)
    x5 = GlobalScope(x5)
    @parameters p1 p2 p3 p4 p5
    p2 = ParentScope(p2)
    p3 = ParentScope(ParentScope(p3))
    p4 = DelayParentScope(p4, 2)
    p5 = GlobalScope(p5)

    j1 = ConstantRateJump(p1, [x1 ~ x1 + 1])
    j2 = MassActionJump(p2, [x2 => 1], [x3 => -1])
    j3 = VariableRateJump(p3, [x3 ~ x3 + 1, x4 ~ x4 + 1])
    j4 = MassActionJump(p4 * p5, [x1 => 1, x5 => 1], [x1 => -1, x5 => -1, x2 => 1])
    @named js = JumpSystem([j1, j2, j3, j4], t, [x1, x2, x3, x4, x5], [p1, p2, p3, p4, p5])

    us = Set()
    ps = Set()
    iv = t
    MT.collect_scoped_vars!(us, ps, js, iv)
    @test issetequal(us, [x2])
    @test issetequal(ps, [p2])

    empty!.((us, ps))
    MT.collect_scoped_vars!(us, ps, js, iv; depth = 0)
    @test issetequal(us, [x1])
    @test issetequal(ps, [p1])

    empty!.((us, ps))
    MT.collect_scoped_vars!(us, ps, js, iv; depth = 1)
    @test issetequal(us, [x2])
    @test issetequal(ps, [p2])

    empty!.((us, ps))
    MT.collect_scoped_vars!(us, ps, js, iv; depth = 2)
    @test issetequal(us, [x3, x4])
    @test issetequal(ps, [p3, p4])

    empty!.((us, ps))
    MT.collect_scoped_vars!(us, ps, js, iv; depth = -1)
    @test issetequal(us, [x5])
    @test issetequal(ps, [p5])
end

# PDMP test
let
    seed = 1111
    Random.seed!(rng, seed)
    @variables X(t) Y(t)
    @parameters k1 k2
    vrj1 = VariableRateJump(k1 * X, [X ~ X - 1]; save_positions = (false, false))
    vrj2 = VariableRateJump(k1, [Y ~ Y + 1]; save_positions = (false, false))
    eqs = [D(X) ~ k2, D(Y) ~ -k2 / 10 * Y]
    @named jsys = JumpSystem([vrj1, vrj2, eqs[1], eqs[2]], t, [X, Y], [k1, k2])
    jsys = complete(jsys)
    X0 = 0.0
    Y0 = 3.0
    u0 = [X => X0, Y => Y0]
    k1val = 1.0
    k2val = 20.0
    p = [k1 => k1val, k2 => k2val]
    tspan = (0.0, 10.0)
    oprob = ODEProblem(jsys, u0, tspan, p)
    jprob = JumpProblem(jsys, oprob; rng, save_positions = (false, false))

    times = range(0.0, tspan[2], length = 100)
    Nsims = 4000
    Xv = zeros(length(times))
    Yv = zeros(length(times))
    for n in 1:Nsims
        sol = solve(jprob, Tsit5(); saveat = times, seed)
        Xv .+= sol[1, :]
        Yv .+= sol[2, :]
        seed += 1
    end
    Xv ./= Nsims
    Yv ./= Nsims

    Xact(t) = X0 * exp(-k1val * t) + (k2val / k1val) * (1 - exp(-k1val * t))
    function Yact(t)
        Y0 * exp(-k2val / 10 * t) + (k1val / (k2val / 10)) * (1 - exp(-k2val / 10 * t))
    end
    @test all(abs.(Xv .- Xact.(times)) .<= 0.05 .* Xv)
    @test all(abs.(Yv .- Yact.(times)) .<= 0.1 .* Yv)
end

# that mixes ODEs and jump types, and then contin events
let
    seed = 1111
    Random.seed!(rng, seed)
    @variables X(t) Y(t)
    @parameters α β
    vrj = VariableRateJump(β * X, [X ~ X - 1]; save_positions = (false, false))
    crj = ConstantRateJump(β * Y, [Y ~ Y - 1])
    maj = MassActionJump(α, [0 => 1], [Y => 1])
    eqs = [D(X) ~ α * (1 + Y)]
    @named jsys = JumpSystem([maj, crj, vrj, eqs[1]], t, [X, Y], [α, β])
    jsys = complete(jsys)
    p = (α = 6.0, β = 2.0, X₀ = 2.0, Y₀ = 1.0)
    u0map = [X => p.X₀, Y => p.Y₀]
    pmap = [α => p.α, β => p.β]
    tspan = (0.0, 20.0)
    oprob = ODEProblem(jsys, u0map, tspan, pmap)
    jprob = JumpProblem(jsys, oprob; rng, save_positions = (false, false))
    times = range(0.0, tspan[2], length = 100)
    Nsims = 4000
    Xv = zeros(length(times))
    Yv = zeros(length(times))
    for n in 1:Nsims
        sol = solve(jprob, Tsit5(); saveat = times, seed)
        Xv .+= sol[1, :]
        Yv .+= sol[2, :]
        seed += 1
    end
    Xv ./= Nsims
    Yv ./= Nsims

    function Yf(t, p)
        local α, β, X₀, Y₀ = p
        return (α / β) + (Y₀ - α / β) * exp(-β * t)
    end
    function Xf(t, p)
        local α, β, X₀, Y₀ = p
        return (α / β) + (α^2 / β^2) + α * (Y₀ - α / β) * t * exp(-β * t) +
               (X₀ - α / β - α^2 / β^2) * exp(-β * t)
    end
    Xact = [Xf(t, p) for t in times]
    Yact = [Yf(t, p) for t in times]
    @test all(abs.(Xv .- Xact) .<= 0.05 .* Xv)
    @test all(abs.(Yv .- Yact) .<= 0.05 .* Yv)

    function affect!(integ, u, p, ctx)
        savevalues!(integ, true)
        terminate!(integ)
        nothing
    end
    cevents = [t ~ 0.2] => (affect!, [], [], [], nothing)
    @named jsys = JumpSystem([maj, crj, vrj, eqs[1]], t, [X, Y], [α, β];
        continuous_events = cevents)
    jsys = complete(jsys)
    tspan = (0.0, 200.0)
    oprob = ODEProblem(jsys, u0map, tspan, pmap)
    jprob = JumpProblem(jsys, oprob; rng, save_positions = (false, false))
    Xsamp = 0.0
    Nsims = 4000
    for n in 1:Nsims
        sol = solve(jprob, Tsit5(); saveat = tspan[2], seed)
        @test sol.retcode == ReturnCode.Terminated
        Xsamp += sol[1, end]
        seed += 1
    end
    Xsamp /= Nsims
    @test abs(Xsamp - Xf(0.2, p) < 0.05 * Xf(0.2, p))
end
