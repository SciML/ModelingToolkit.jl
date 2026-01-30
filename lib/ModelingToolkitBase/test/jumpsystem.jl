using ModelingToolkitBase, DiffEqBase, JumpProcesses, Test, LinearAlgebra
using SymbolicIndexingInterface, OrderedCollections
using Random, StableRNGs, NonlinearSolve
using OrdinaryDiffEq
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using BenchmarkTools
using Symbolics: SymbolicT, unwrap
import SymbolicIndexingInterface as SII
MT = ModelingToolkitBase

rng = StableRNG(12345)

# basic MT SIR model with tweaks
@parameters β γ
@constants h = 1
@variables S(t) I(t) R(t)
rate₁ = β * S * I * h
affect₁ = [S ~ Pre(S) - 1 * h, I ~ Pre(I) + 1]
rate₂ = γ * I + t
affect₂ = [I ~ Pre(I) - 1, R ~ Pre(R) + 1]
j₁ = ConstantRateJump(rate₁, affect₁)
j₂ = VariableRateJump(rate₂, affect₂)
@named js = JumpSystem([j₁, j₂], t, [S, I, R], [β, γ, h])
unknowntoid = Dict(MT.value(unknown) => i for (i, unknown) in enumerate(unknowns(js)))
mtjump1 = MT.assemble_crj(js, j₁, unknowntoid)
mtjump2 = MT.assemble_vrj(js, j₂, unknowntoid)

# doc version
rate1(u, p, t) = (0.1 / 1000.0) * u[1] * u[2]
function affect1!(integrator)
    integrator.u[1] -= 1
    return integrator.u[2] += 1
end
jump1 = ConstantRateJump(rate1, affect1!)
rate2(u, p, t) = 0.01u[2] + t
function affect2!(integrator)
    integrator.u[2] -= 1
    return integrator.u[3] += 1
end
jump2 = VariableRateJump(rate2, affect2!)

# test crjs
u = [100, 9, 5]
p = (0.1 / 1000, 0.01, 1)
tf = 1.0
mutable struct TestInt{U, V, T}
    u::U
    p::V
    t::T
end
SII.state_values(x::TestInt) = x.u
SII.parameter_values(x::TestInt) = x.p
SII.current_time(x::TestInt) = x.t
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
affect₃ = [I ~ Pre(I) * h - 1, R ~ Pre(R) + 1]
j₃ = ConstantRateJump(rate₃, affect₃)
@named js2 = JumpSystem([j₁, j₃], t, [S, I, R], [β, γ, h])
js2 = complete(js2)
u₀ = [999, 1, 0];
tspan = (0.0, 250.0);
u₀map = [S => 999, I => 1, R => 0]
parammap = [β => 0.1 / 1000, γ => 0.01]
jprob = JumpProblem(
    js2, [u₀map; parammap], tspan; aggregator = Direct(),
    save_positions = (false, false), rng
)
p = parameter_values(jprob)
@test jprob.prob isa DiscreteProblem
Nsims = 1000
function getmean(jprob, Nsims; use_stepper = true)
    m = 0.0
    for i in 1:Nsims
        i % 200 == 0 && @info i
        sol = use_stepper ? solve(jprob, SSAStepper()) : solve(jprob)
        m += sol[end, end]
    end
    return m / Nsims
end
m = getmean(jprob, Nsims)

# test auto-alg selection works
jprobb = JumpProblem(js2, [u₀map; parammap], tspan; save_positions = (false, false), rng)
mb = getmean(jprobb, Nsims; use_stepper = false)
@test abs(m - mb) / m < 0.01

@variables S2(t)
obs = [S2 ~ 2 * S]
@named js2b = JumpSystem([j₁, j₃], t, [S, I, R], [β, γ, h], observed = obs)
js2b = complete(js2b)
jprob = JumpProblem(
    js2b, [u₀map; parammap], tspan; aggregator = Direct(),
    save_positions = (false, false), rng
)
@test jprob.prob isa DiscreteProblem
sol = solve(jprob, SSAStepper(); saveat = tspan[2] / 10)
@test all(2 .* sol[S] .== sol[S2])

# test save_positions is working
jprob = JumpProblem(
    js2, [u₀map; parammap], tspan; aggregator = Direct(),
    save_positions = (false, false), rng
)
sol = solve(jprob, SSAStepper(); saveat = 1.0)
@test all((sol.t) .== collect(0.0:tspan[2]))

#test the MT JumpProblem rates/affects are correct
rate2(u, p, t) = 0.01u[2]
jump2 = ConstantRateJump(rate2, affect2!)
mtjumps = jprob.discrete_jump_aggregation
@test abs(mtjumps.rates[1](u, p, tf) - jump1.rate(u, p, tf)) < 10 * eps()
@test abs(mtjumps.rates[2](u, p, tf) - jump2.rate(u, p, tf)) < 10 * eps()

ModelingToolkitBase.@set! mtintegrator.p = parameter_values(jprob)
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
    return integrator.u[2] += 1
end
j1 = ConstantRateJump(r1, a1!)
r2(u, p, t) = 0.01u[2]
function a2!(integrator)
    integrator.u[2] -= 1
    return integrator.u[3] += 1
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
jprob = JumpProblem(js3, [u₀map; parammap], tspan; aggregator = Direct(), rng)
@test jprob.prob isa DiscreteProblem
m3 = getmean(jprob, Nsims)
@test abs(m - m3) / m < 0.01

# maj jump test with various dep graphs
@named js3b = JumpSystem([maj1, maj2], t, [S, I, R], [β, γ])
js3b = complete(js3b)
jprobb = JumpProblem(js3b, [u₀map; parammap], tspan; aggregator = NRM(), rng)
@test jprobb.prob isa DiscreteProblem
m4 = getmean(jprobb, Nsims)
@test abs(m - m4) / m < 0.01
jprobc = JumpProblem(js3b, [u₀map; parammap], tspan; aggregator = RSSA(), rng)
@test jprobc.prob isa DiscreteProblem
m4 = getmean(jprobc, Nsims)
@test abs(m - m4) / m < 0.01

# mass action jump tests for other reaction types (zero order, decay)
maj1 = MassActionJump(2.0, [0 => 1], [S => 1])
maj2 = MassActionJump(γ, [S => 1], [S => -1])
@named js4 = JumpSystem([maj1, maj2], t, [S], [β, γ])
js4 = complete(js4)
jprob = JumpProblem(
    js4, [S => 999, β => 100.0, γ => 0.01], (0, 1000.0); aggregator = Direct(), rng
)
@test jprob.prob isa DiscreteProblem
m4 = getmean(jprob, Nsims)
@test abs(m4 - 2.0 / 0.01) * 0.01 / 2.0 < 0.01

# test second order rx runs
maj1 = MassActionJump(2.0, [0 => 1], [S => 1])
maj2 = MassActionJump(γ, [S => 2], [S => -1])
@named js4 = JumpSystem([maj1, maj2], t, [S], [β, γ])
js4 = complete(js4)
jprob = JumpProblem(
    js4, [S => 999, β => 100.0, γ => 0.01], (0, 1000.0); aggregator = Direct(), rng
)
@test jprob.prob isa DiscreteProblem
sol = solve(jprob, SSAStepper());

# issue #819
@testset "Combined system name collisions" begin
    sys1 = JumpSystem([maj1, maj2], t, [S], [β, γ], name = :sys1)
    sys2 = JumpSystem([maj1, maj2], t, [S], [β, γ], name = :sys1)
    @test_throws ModelingToolkitBase.NonUniqueSubsystemsError JumpSystem(
        [sys1.γ ~ sys2.γ], t, [], [],
        systems = [sys1, sys2], name = :foo
    )
end

# test if param mapper is setup correctly for callbacks
@testset "Parammapper with callbacks" begin
    @parameters k1 k2 k3
    @variables A(t) B(t)
    maj1 = MassActionJump(k1 * k3, [0 => 1], [A => -1, B => 1])
    maj2 = MassActionJump(k2, [B => 1], [A => 1, B => -1])
    @named js5 = JumpSystem([maj1, maj2], t, [A, B], [k1, k2, k3])
    js5 = complete(js5)
    p = [k1 => 2.0, k2 => 0.0, k3 => 0.5]
    u₀ = [A => 100, B => 0]
    tspan = (0.0, 2000.0)
    jprob = JumpProblem(
        js5, [u₀; p], tspan; aggregator = Direct(), save_positions = (false, false), rng
    )
    @test jprob.prob isa DiscreteProblem
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
@testset "Observed handling tests" begin
    @variables OBS(t)
    @named js5 = JumpSystem([maj1, maj2], t, [S], [β, γ]; observed = [OBS ~ 2 * S * h])
    OBS2 = OBS
    @test isequal(OBS2, @nonamespace js5.OBS)
    @unpack OBS = js5
    @test isequal(OBS2, OBS)
end

# test to make sure dep graphs are correct
@testset "Dependency graph tests" begin
    # A + 2X --> 3X
    # 3X --> A + 2X
    # B --> X
    # X --> B
    @variables A(t) X(t) B(t)
    jumps = [
        MassActionJump(1.0, [A => 1, X => 2], [A => -1, X => 1]),
        MassActionJump(1.0, [X => 3], [A => 1, X => -1]),
        MassActionJump(1.0, [B => 1], [B => -1, X => 1]),
        MassActionJump(1.0, [X => 1], [B => 1, X => -1]),
    ]
    @named js = JumpSystem(jumps, t, [A, X, B], [])
    jdeps = asgraph(js; eqs = MT.jumps(js))
    vdeps = variable_dependencies(js; eqs = MT.jumps(js))
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

crj = ConstantRateJump(1.0, [X ~ Pre(X) - 1])
js1 = complete(JumpSystem([crj], t, [X], [k]; name = :js1))
js2 = complete(JumpSystem([crj], t, [X], []; name = :js2))

maj = MassActionJump(1.0, [X => 1], [X => -1])
js3 = complete(JumpSystem([maj], t, [X], [k]; name = :js2))
js4 = complete(JumpSystem([maj], t, [X], []; name = :js3))

u0 = [X => 10]
tspan = (0.0, 1.0)
ps = [k => 1.0]

@test_nowarn jp1 = JumpProblem(js1, [u0; ps], tspan; aggregator = Direct())
@test_nowarn jp2 = JumpProblem(js2, u0, tspan; aggregator = Direct())
@test_nowarn jp3 = JumpProblem(js3, [u0; ps], tspan; aggregator = Direct())
@test_nowarn jp4 = JumpProblem(js4, u0, tspan; aggregator = Direct())

# Ensure `mtkcompile` (and `@mtkcompile`) works on JumpSystem (by doing nothing)
# Issue#2558
@parameters k
@variables X(t)
rate = k
affect = [X ~ Pre(X) - 1]

j1 = ConstantRateJump(k, [X ~ Pre(X) - 1])
@test_nowarn @mtkcompile js1 = JumpSystem([j1], t, [X], [k])

# test correct autosolver is selected, which implies appropriate dep graphs are available
@testset "Autosolver test" begin
    @parameters k
    @variables X(t)
    rate = k
    affect = [X ~ Pre(X) - 1]
    j1 = ConstantRateJump(k, [X ~ Pre(X) - 1])

    Nv = [1, JumpProcesses.USE_DIRECT_THRESHOLD + 1, JumpProcesses.USE_RSSA_THRESHOLD + 1]
    algtypes = [Direct, RSSA, RSSACR]
    for (N, algtype) in zip(Nv, algtypes)
        @named jsys = JumpSystem([deepcopy(j1) for _ in 1:N], t, [X], [k])
        jsys = complete(jsys)
        jprob = JumpProblem(jsys, [X => 10, k => 1], (0.0, 10.0))
        @test jprob.aggregator isa algtype
    end
end

# basic VariableRateJump test
@testset "VRJ test" begin
    N = 1000  # number of simulations for testing solve accuracy
    Random.seed!(rng, 1111)
    @variables A(t) B(t) C(t)
    @parameters k
    vrj = VariableRateJump(k * (sin(t) + 1), [A ~ Pre(A) + 1, C ~ Pre(C) + 2])
    js = complete(JumpSystem([vrj], t, [A, C], [k]; name = :js, observed = [B ~ C * A]))
    jprob = JumpProblem(
        js, [A => 0, C => 0, k => 1], (0.0, 10.0); aggregator = Direct(), rng
    )
    @test jprob.prob isa ODEProblem
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
@testset "`collect_vars!` for jumps" begin
    @variables x1(t) x2(t) x3(t) x4(t) x5(t)
    @parameters p1 p2 p3 p4 p5
    j1 = ConstantRateJump(p1, [x1 ~ Pre(x1) + 1])
    j2 = MassActionJump(p2, [x2 => 1], [x3 => -1])
    j3 = VariableRateJump(p3, [x3 ~ Pre(x3) + 1, x4 ~ Pre(x4) + 1])
    j4 = MassActionJump(p4 * p5, [x1 => 1, x5 => 1], [x1 => -1, x5 => -1, x2 => 1])
    us = OrderedSet{SymbolicT}()
    ps = OrderedSet{SymbolicT}()
    iv = unwrap(t)

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
@testset "Scoping tests" begin
    @variables x1(t) x2(t) x3(t) x4(t)
    x2 = ParentScope(x2)
    x3 = ParentScope(ParentScope(x3))
    x4 = GlobalScope(x4)
    @parameters p1 p2 p3 p4
    p2 = ParentScope(p2)
    p3 = ParentScope(ParentScope(p3))
    p4 = GlobalScope(p4)

    j1 = ConstantRateJump(p1, [x1 ~ Pre(x1) + 1])
    j2 = MassActionJump(p2, [x2 => 1], [x3 => -1])
    j3 = VariableRateJump(p3, [x3 ~ Pre(x3) + 1, x4 ~ Pre(x4) + 1])
    j4 = MassActionJump(p4 * p4, [x1 => 1, x4 => 1], [x1 => -1, x4 => -1, x2 => 1])
    @named js = JumpSystem([j1, j2, j3, j4], t, [x1, x2, x3, x4], [p1, p2, p3, p4])

    us = OrderedSet{SymbolicT}()
    ps = OrderedSet{SymbolicT}()
    iv = unwrap(t)
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
    @test issetequal(us, [x3])
    @test issetequal(ps, [p3])

    empty!.((us, ps))
    MT.collect_scoped_vars!(us, ps, js, iv; depth = -1)
    @test issetequal(us, [x4])
    @test issetequal(ps, [p4])
end

# PDMP test
@testset "PDMP test" begin
    seed = 1111
    Random.seed!(rng, seed)
    @variables X(t) Y(t)
    @parameters k1 k2
    vrj1 = VariableRateJump(k1 * X, [X ~ Pre(X) - 1]; save_positions = (false, false))
    vrj2 = VariableRateJump(k1, [Y ~ Pre(Y) + 1]; save_positions = (false, false))
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
    jprob = JumpProblem(jsys, [u0; p], tspan; rng, save_positions = (false, false))

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
@testset "ODEs + Jumps + Continuous events" begin
    seed = 1111
    Random.seed!(rng, seed)
    @variables X(t) Y(t)
    @parameters α β
    vrj = VariableRateJump(β * X, [X ~ Pre(X) - 1]; save_positions = (false, false))
    crj = ConstantRateJump(β * Y, [Y ~ Pre(Y) - 1])
    maj = MassActionJump(α, [0 => 1], [Y => 1])
    eqs = [D(X) ~ α * (1 + Y)]
    @named jsys = JumpSystem([maj, crj, vrj, eqs[1]], t, [X, Y], [α, β])
    jsys = complete(jsys)
    p = (α = 6.0, β = 2.0, X₀ = 2.0, Y₀ = 1.0)
    u0map = [X => p.X₀, Y => p.Y₀]
    pmap = [α => p.α, β => p.β]
    tspan = (0.0, 20.0)
    jprob = JumpProblem(jsys, [u0map; pmap], tspan; rng, save_positions = (false, false))
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

    function affect!(mod, obs, ctx, integ)
        savevalues!(integ, true)
        terminate!(integ)
        (;)
    end
    cevents = [t ~ 0.2] => (; f = affect!)
    @named jsys = JumpSystem(
        [maj, crj, vrj, eqs[1]], t, [X, Y], [α, β];
        continuous_events = cevents
    )
    jsys = complete(jsys)
    tspan = (0.0, 200.0)
    jprob = JumpProblem(jsys, [u0map; pmap], tspan; rng, save_positions = (false, false))
    Xsamp = 0.0
    Nsims = 4000
    for n in 1:Nsims
        sol = solve(jprob, Tsit5(); saveat = tspan[2], seed)
        @test sol.retcode == ReturnCode.Terminated
        Xsamp += sol[1, end]
        seed += 1
    end
    Xsamp /= Nsims
    @test abs(Xsamp - Xf(0.2, p)) < 0.05 * Xf(0.2, p)
end

@testset "JumpProcess simulation should be Int64 valued (#3446)" begin
    @parameters p d
    @variables X(t)
    rate1 = p
    rate2 = X * d
    affect1 = [X ~ Pre(X) + 1]
    affect2 = [X ~ Pre(X) - 1]
    j1 = ConstantRateJump(rate1, affect1)
    j2 = ConstantRateJump(rate2, affect2)

    # Works.
    @mtkcompile js = JumpSystem([j1, j2], t, [X], [p, d])
    jprob = JumpProblem(
        js, [X => 15, p => 2.0, d => 0.5], (0.0, 10.0); aggregator = Direct(), u0_eltype = Int
    )
    sol = solve(jprob, SSAStepper())
    @test eltype(sol[X]) === Int64
end

@testset "Issue#3571: `remake(::JumpProblem)`" begin
    @variables X(t)
    @parameters a b
    eq = D(X) ~ a
    rate = b * X
    affect = [X ~ Pre(X) - 1]
    crj = ConstantRateJump(rate, affect)
    @named jsys = JumpSystem([crj, eq], t, [X], [a, b])
    jsys = complete(jsys)
    jprob = JumpProblem(jsys, [:X => 1.0, :a => 1.0, :b => 0.5], (0.0, 10.0))
    jprob2 = remake(jprob; u0 = [:X => 10.0])
    @test jprob2[X] ≈ 10.0
end

@testset "Proper substitution in `JumpSysMajParamWrapper`" begin
    @variables X(t)
    @parameters p d
    jump1 = MassActionJump(exp(p), Pair{Num, Real}[], [X => 1], nothing)
    jump2 = ConstantRateJump(d * exp(X) * X, [X ~ Pre(X) - 1])
    @named sys = JumpSystem([jump1, jump2], t, [X], [p, d])
    sys = complete(sys)
    @test_nowarn JumpProblem(sys, [:X => 0.1, :p => 1.0, :d => 2.0], (0.0, 1.0))
end

@testset "Issue#4216: Discrete events fire once in hybrid JumpProblems" begin
    @variables X(t)
    @parameters a b

    # Create ODE equation
    eq = D(X) ~ a

    # Create jump
    crj = ConstantRateJump(b * X, [X ~ Pre(X) - 1])

    # Discrete event: at t=1.0, add 5.0 to X
    discrete_event = [1.0] => [X ~ Pre(X) + 5.0]

    # Create hybrid JumpSystem with ODE + jump + discrete event
    @named jsys = JumpSystem([crj, eq], t, [X], [a, b]; discrete_events = [discrete_event])
    jsys = complete(jsys)

    # Create JumpProblem
    jprob = JumpProblem(jsys, [:X => 10.0, :a => 1.0, :b => 0.01], (0.0, 3.0); rng)

    # Callback should NOT be stored in the underlying problem (only in JumpProblem)
    @test !haskey(jprob.prob.kwargs, :callback)
    @test haskey(jprob.kwargs, :callback)

    # Solve and verify discrete event fires exactly once
    sol = solve(jprob, Tsit5())
    X_before = sol(0.99; idxs = X)
    X_after = sol(1.01; idxs = X)
    change = X_after - X_before

    # The change should be approximately 5.0 (not ~10.0 which would indicate double firing)
    @test isapprox(change, 5.0, atol = 0.1)
end

@testset "Issue#4216: Continuous events not duplicated in hybrid JumpProblems" begin
    @variables X(t)
    @parameters a b

    # Create ODE equation: X grows linearly
    eq = D(X) ~ a

    # Create jump
    crj = ConstantRateJump(b * X, [X ~ Pre(X) - 1])

    # Continuous event: when X crosses 15.0, add 5.0 to X
    continuous_event = [X ~ 15.0] => [X ~ Pre(X) + 5.0]

    # Create hybrid JumpSystem with ODE + jump + continuous event
    @named jsys = JumpSystem(
        [crj, eq], t, [X], [a, b]; continuous_events = [continuous_event]
    )
    jsys = complete(jsys)

    # Create JumpProblem (starting at X=10, with a=10 so X grows quickly)
    jprob = JumpProblem(jsys, [:X => 10.0, :a => 10.0, :b => 0.001], (0.0, 3.0); rng)

    # Callback should NOT be stored in the underlying problem (only in JumpProblem)
    @test !haskey(jprob.prob.kwargs, :callback)
    @test haskey(jprob.kwargs, :callback)

    # Solve - the continuous event should fire when X crosses 15.0
    sol = solve(jprob, Tsit5())

    # After crossing 15.0, X should jump by 5.0 (to ~20.0), not by 10.0 (which would be ~25.0)
    # With a=10, X reaches 15 at t≈0.5, then jumps to ~20, continues growing
    # At t=1.0, X should be around 20 + 5 = 25 (not 30 which would indicate double firing)
    X_at_1 = sol(1.0; idxs = X)

    # X starts at 10, grows at rate 10, crosses 15 at t≈0.5, jumps by 5 to ~20
    # Then continues growing: at t=1.0, X ≈ 20 + 10*(1.0-0.5) = 25
    # If event fired twice: X ≈ 25 + 10*(1.0-0.5) = 30
    @test X_at_1 < 28.0  # Should be ~25, not ~30
end

# Test save_positions kwarg is correctly forwarded to JumpProcesses for discrete jumps
# Note: save_positions to JumpProblem controls MAJs and CRJs only.
# VRJs have their own save_positions set at construction time.
@testset "save_positions kwarg forwarding" begin
    # Test 1: DiscreteProblem-based JumpProblem with MassActionJumps
    @testset "MassActionJump with DiscreteProblem" begin
        @variables A(t)
        @parameters k
        maj = MassActionJump(k, [A => 1], [A => -1])
        @named jsys = JumpSystem([maj], t, [A], [k])
        jsys = complete(jsys)

        # With save_positions=(false, false) and saveat, should get exact number of points
        jprob = JumpProblem(
            jsys, [A => 100, k => 1.0], (0.0, 10.0);
            aggregator = Direct(), save_positions = (false, false), rng
        )
        @test jprob.prob isa DiscreteProblem

        # Verify save_positions reaches the aggregator
        @test jprob.discrete_jump_aggregation.save_positions == (false, false)

        # Solve with saveat and verify no extra points from jumps
        times = 0.0:1.0:10.0
        sol = solve(jprob, SSAStepper(); saveat = times)
        @test length(sol.t) == length(times)
        @test all(sol.t .== collect(times))
    end

    # Test 2: DiscreteProblem-based JumpProblem with ConstantRateJumps
    @testset "ConstantRateJump with DiscreteProblem" begin
        @variables A(t)
        @parameters k
        crj = ConstantRateJump(k * A, [A ~ Pre(A) - 1])
        @named jsys = JumpSystem([crj], t, [A], [k])
        jsys = complete(jsys)

        jprob = JumpProblem(
            jsys, [A => 100, k => 0.1], (0.0, 10.0);
            aggregator = Direct(), save_positions = (false, false), rng
        )
        @test jprob.prob isa DiscreteProblem

        # Verify save_positions reaches the aggregator
        @test jprob.discrete_jump_aggregation.save_positions == (false, false)

        times = 0.0:1.0:10.0
        sol = solve(jprob, SSAStepper(); saveat = times)
        @test length(sol.t) == length(times)
    end

    # Test 3: ODEProblem-based JumpProblem with VariableRateJumps only
    # VRJs have their own save_positions - the JumpProblem-level save_positions
    # should not cause an error (i.e., should not be passed to ODEProblem)
    @testset "VariableRateJump with ODEProblem" begin
        @variables A(t)
        @parameters k
        vrj = VariableRateJump(k * (1 + sin(t)), [A ~ Pre(A) + 1])
        @named jsys = JumpSystem([vrj], t, [A], [k])
        jsys = complete(jsys)

        # This previously errored because save_positions was passed to ODEProblem
        jprob = JumpProblem(
            jsys, [A => 0, k => 1.0], (0.0, 10.0);
            aggregator = Direct(), save_positions = (false, false), rng
        )
        @test jprob.prob isa ODEProblem

        # Verify save_positions is NOT in the ODEProblem kwargs (this was the bug)
        @test !haskey(jprob.prob.kwargs, :save_positions)

        # Should solve without error
        sol = solve(jprob, Tsit5())
        @test SciMLBase.successful_retcode(sol)
    end

    # Test 4: Hybrid system with ODEs and ConstantRateJumps
    @testset "Hybrid ODE + ConstantRateJump" begin
        @variables X(t)
        @parameters a b
        eq = D(X) ~ a
        crj = ConstantRateJump(b * X, [X ~ Pre(X) - 1])
        @named jsys = JumpSystem([crj, eq], t, [X], [a, b])
        jsys = complete(jsys)

        jprob = JumpProblem(
            jsys, [X => 10.0, a => 1.0, b => 0.1], (0.0, 10.0);
            save_positions = (false, false), rng
        )
        @test jprob.prob isa ODEProblem

        # Verify save_positions is NOT in the ODEProblem kwargs
        @test !haskey(jprob.prob.kwargs, :save_positions)

        # Verify save_positions reaches the discrete aggregator
        @test jprob.discrete_jump_aggregation.save_positions == (false, false)

        sol = solve(jprob, Tsit5())
        @test SciMLBase.successful_retcode(sol)

        times = 0.0:1.0:10.0
        sol = solve(jprob, Tsit5(); saveat = times)
        @test length(sol.t) == length(times)
    end

    # Test 5: Hybrid system with ODEs, VRJs, CRJs, and MAJs
    # save_positions should reach the discrete aggregator for CRJs/MAJs
    # VRJs use their own save_positions from construction
    @testset "Hybrid ODE + all jump types" begin
        @variables X(t) Y(t)
        @parameters a b c d
        eq = D(X) ~ a
        # VRJ with its own save_positions
        vrj = VariableRateJump(b * (1 + sin(t)), [X ~ Pre(X) + 1]; save_positions = (false, false))
        crj = ConstantRateJump(c * Y, [Y ~ Pre(Y) - 1])
        maj = MassActionJump(d, [0 => 1], [Y => 1])
        @named jsys = JumpSystem([vrj, crj, maj, eq], t, [X, Y], [a, b, c, d])
        jsys = complete(jsys)

        jprob = JumpProblem(
            jsys, [X => 0.0, Y => 10, a => 0.1, b => 0.5, c => 0.1, d => 1.0], (0.0, 10.0);
            save_positions = (false, false), rng
        )
        @test jprob.prob isa ODEProblem

        # Verify save_positions is NOT in the ODEProblem kwargs
        @test !haskey(jprob.prob.kwargs, :save_positions)

        # Verify save_positions reaches the discrete aggregator (for CRJs/MAJs)
        @test jprob.discrete_jump_aggregation.save_positions == (false, false)

        sol = solve(jprob, Tsit5())
        @test SciMLBase.successful_retcode(sol)

        # With all jumps having save_positions=(false,false), saveat should give exact points
        times = 0.0:1.0:10.0
        sol = solve(jprob, Tsit5(); saveat = times)
        @test length(sol.t) == length(times)
    end

    # Test 6: Default save_positions for discrete jumps should be (true, true)
    @testset "Default save_positions for discrete jumps" begin
        @variables A(t)
        @parameters k
        maj = MassActionJump(k, [A => 1], [A => -1])
        @named jsys = JumpSystem([maj], t, [A], [k])
        jsys = complete(jsys)

        # No save_positions specified - should default to (true, true) for discrete aggregator
        jprob = JumpProblem(jsys, [A => 100, k => 1.0], (0.0, 10.0); rng)

        @test jprob.discrete_jump_aggregation.save_positions == (true, true)
    end
end

# Test Symbol keys work with discrete events in JumpProblems
# This tests the fix for the bug where Symbol keys in op failed to convert
# when process_events created ImplicitDiscreteProblems for affect subsystems
@testset "Symbol keys with discrete events in JumpProblem" begin
    # Test 1: Pure discrete jump system with Symbol keys and discrete event
    # The main test is that JumpProblem construction doesn't error with Symbol keys
    @testset "Pure discrete jumps with Symbol keys and discrete events" begin
        @variables X(t)
        @parameters k X0

        # A discrete event that resets X to parameter X0
        discrete_event = [5.0] => [X ~ X0]

        # Mass action jump: X decays
        maj = MassActionJump(k, [X => 1], [X => -1])

        @named jsys = JumpSystem([maj], t, [X], [k, X0]; discrete_events = [discrete_event])
        jsys = complete(jsys)

        # Using Symbol keys - this previously failed with TypeError because :X0
        # couldn't be found in the affect system's symbol table
        jprob = JumpProblem(
            jsys, [:X => 50, :k => 0.1, :X0 => 100], (0.0, 10.0);
            aggregator = Direct(), rng
        )
        @test jprob.prob isa DiscreteProblem

        # Solve and verify the discrete event fires correctly
        sol = solve(jprob, SSAStepper())
        @test SciMLBase.successful_retcode(sol)

        # After t=5.0, X should have been reset to X0=100
        # (may be off by a few due to jumps occurring after the reset)
        idx_after = findfirst(t -> t > 5.0, sol.t)
        if idx_after !== nothing
            # X should be close to 100 (reset value), allow for a few jumps
            @test sol[X][idx_after] >= 95  # reset to 100, minus a few jumps
        end
    end

    # Test 2: Constant rate jump with Symbol keys and discrete event
    @testset "ConstantRateJump with Symbol keys and discrete events" begin
        @variables X(t)
        @parameters k reset_val

        discrete_event = [2.0] => [X ~ reset_val]
        crj = ConstantRateJump(k * X, [X ~ Pre(X) - 1])

        @named jsys = JumpSystem([crj], t, [X], [k, reset_val]; discrete_events = [discrete_event])
        jsys = complete(jsys)

        # Using Symbol keys - main test is that this doesn't error
        jprob = JumpProblem(
            jsys, [:X => 20, :k => 0.5, :reset_val => 50], (0.0, 5.0);
            aggregator = Direct(), rng
        )
        @test jprob.prob isa DiscreteProblem

        sol = solve(jprob, SSAStepper())
        @test SciMLBase.successful_retcode(sol)

        # Verify X was reset at t=2.0 (allow for a few jumps)
        idx_after = findfirst(t -> t > 2.0, sol.t)
        if idx_after !== nothing
            @test sol[X][idx_after] >= 45  # reset to 50, minus a few jumps
        end
    end

    # Test 3: Hybrid system (ODE + jumps) with Symbol keys and discrete events
    @testset "Hybrid ODE + jumps with Symbol keys and discrete events" begin
        @variables X(t)
        @parameters a b X0

        eq = D(X) ~ a
        crj = ConstantRateJump(b * X, [X ~ Pre(X) - 1])
        discrete_event = [1.5] => [X ~ X0]

        @named jsys = JumpSystem([crj, eq], t, [X], [a, b, X0]; discrete_events = [discrete_event])
        jsys = complete(jsys)

        # Using Symbol keys
        jprob = JumpProblem(
            jsys, [:X => 10.0, :a => 1.0, :b => 0.01, :X0 => 50.0], (0.0, 3.0); rng
        )
        @test jprob.prob isa ODEProblem

        sol = solve(jprob, Tsit5())
        @test SciMLBase.successful_retcode(sol)

        # Verify the reset happened - X should be ~50 just after t=1.5
        X_after = sol(1.51; idxs = X)
        @test isapprox(X_after, 50.0, atol = 1.0)
    end

    # Test 4: Multiple discrete events with Symbol keys
    @testset "Multiple discrete events with Symbol keys" begin
        @variables X(t)
        @parameters k val1 val2

        event1 = [2.0] => [X ~ val1]
        event2 = [4.0] => [X ~ val2]
        maj = MassActionJump(k, [X => 1], [X => -1])

        @named jsys = JumpSystem([maj], t, [X], [k, val1, val2]; discrete_events = [event1, event2])
        jsys = complete(jsys)

        # Main test: construction with Symbol keys doesn't error
        jprob = JumpProblem(
            jsys, [:X => 100, :k => 0.05, :val1 => 200, :val2 => 50], (0.0, 6.0);
            aggregator = Direct(), rng
        )

        sol = solve(jprob, SSAStepper())
        @test SciMLBase.successful_retcode(sol)

        # Check first reset at t=2 (allow for jumps to occur)
        idx1 = findfirst(t -> t > 2.0, sol.t)
        if idx1 !== nothing
            @test sol[X][idx1] >= 195  # reset to 200, minus a few jumps
        end

        # Check second reset at t=4 (allow for jumps)
        idx2 = findfirst(t -> t > 4.0, sol.t)
        if idx2 !== nothing
            @test sol[X][idx2] >= 45  # reset to 50, minus a few jumps
        end
    end

    # Test 5: VariableRateJump with Symbol keys and discrete events
    @testset "VariableRateJump with Symbol keys and discrete events" begin
        @variables X(t)
        @parameters k X0

        vrj = VariableRateJump(k * (1 + sin(t)), [X ~ Pre(X) + 1])
        discrete_event = [3.0] => [X ~ X0]

        @named jsys = JumpSystem([vrj], t, [X], [k, X0]; discrete_events = [discrete_event])
        jsys = complete(jsys)

        # VRJ systems get an ODEProblem
        jprob = JumpProblem(
            jsys, [:X => 0, :k => 1.0, :X0 => 100], (0.0, 5.0); rng
        )
        @test jprob.prob isa ODEProblem

        sol = solve(jprob, Tsit5())
        @test SciMLBase.successful_retcode(sol)

        # Verify reset at t=3
        X_after = sol(3.01; idxs = X)
        @test isapprox(X_after, 100.0, atol = 2.0)
    end

    # Test 6: Continuous events with Symbol keys
    # The fix also affects continuous events since process_events handles both types
    @testset "Continuous events with Symbol keys" begin
        @variables X(t)
        @parameters a b threshold

        eq = D(X) ~ a
        crj = ConstantRateJump(b * X, [X ~ Pre(X) - 1])

        # Continuous event with parameter in condition
        continuous_event = [X ~ threshold] => [X ~ Pre(X) + 5.0]

        @named jsys = JumpSystem(
            [crj, eq], t, [X], [a, b, threshold];
            continuous_events = [continuous_event]
        )
        jsys = complete(jsys)

        # Using Symbol keys - main test is that JumpProblem construction doesn't error
        # Previously this would fail with TypeError when Symbol keys couldn't be
        # converted in compile_equational_affect
        jprob = JumpProblem(
            jsys, [:X => 10.0, :a => 10.0, :b => 0.001, :threshold => 15.0],
            (0.0, 3.0); rng
        )
        @test jprob.prob isa ODEProblem

        # Callbacks should be in the JumpProblem (same pattern as existing tests)
        @test !haskey(jprob.prob.kwargs, :callback)
        @test haskey(jprob.kwargs, :callback)

        sol = solve(jprob, Tsit5())
        @test SciMLBase.successful_retcode(sol)

        # Correctness test matching Issue#4216 pattern:
        # X starts at 10, grows at rate 10, crosses 15 at t≈0.5
        # If event fires once: X jumps to ~20, continues to ~25 at t=1
        # If event fires twice: X would be ~30 at t=1
        # If event doesn't fire: X = 10 + 10*1 = 20 at t=1
        X_at_1 = sol(1.0; idxs = X)
        @test X_at_1 < 28.0  # Should not fire twice (same check as Issue#4216)
    end
end
