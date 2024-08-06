using ModelingToolkit, DiffEqBase, JumpProcesses, Test, LinearAlgebra, StableRNGs
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
jprob = JumpProblem(js2, dprob, Direct(), save_positions = (false, false), rng = rng)
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
jprob = JumpProblem(js2b, dprob, Direct(), save_positions = (false, false), rng = rng)
sol = solve(jprob, SSAStepper(), saveat = tspan[2] / 10)
@test all(2 .* sol[S] .== sol[S2])

# test save_positions is working
jprob = JumpProblem(js2, dprob, Direct(), save_positions = (false, false), rng = rng)
sol = solve(jprob, SSAStepper(), saveat = 1.0)
@test all((sol.t) .== collect(0.0:tspan[2]))

#test the MT JumpProblem rates/affects are correct
rate2(u, p, t) = 0.01u[2]
jump2 = ConstantRateJump(rate2, affect2!)
mtjumps = jprob.discrete_jump_aggregation
@test abs(mtjumps.rates[1](u, p, tf) - jump1.rate(u, p, tf)) < 10 * eps()
@test abs(mtjumps.rates[2](u, p, tf) - jump2.rate(u, p, tf)) < 10 * eps()
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
jprob = JumpProblem(prob, Direct(), jset, save_positions = (false, false), rng = rng)
m2 = getmean(jprob, Nsims)

# test JumpSystem solution agrees with direct version
@test abs(m - m2) / m < 0.01

# mass action jump tests for SIR model
maj1 = MassActionJump(2 * β / 2, [S => 1, I => 1], [S => -1, I => 1])
maj2 = MassActionJump(γ, [I => 1], [I => -1, R => 1])
@named js3 = JumpSystem([maj1, maj2], t, [S, I, R], [β, γ])
js3 = complete(js3)
dprob = DiscreteProblem(js3, u₀map, tspan, parammap)
jprob = JumpProblem(js3, dprob, Direct(), rng = rng)
m3 = getmean(jprob, Nsims)
@test abs(m - m3) / m < 0.01

# maj jump test with various dep graphs
@named js3b = JumpSystem([maj1, maj2], t, [S, I, R], [β, γ])
js3b = complete(js3b)
jprobb = JumpProblem(js3b, dprob, NRM(), rng = rng)
m4 = getmean(jprobb, Nsims)
@test abs(m - m4) / m < 0.01
jprobc = JumpProblem(js3b, dprob, RSSA(), rng = rng)
m4 = getmean(jprobc, Nsims)
@test abs(m - m4) / m < 0.01

# mass action jump tests for other reaction types (zero order, decay)
maj1 = MassActionJump(2.0, [0 => 1], [S => 1])
maj2 = MassActionJump(γ, [S => 1], [S => -1])
@named js4 = JumpSystem([maj1, maj2], t, [S], [β, γ])
js4 = complete(js4)
dprob = DiscreteProblem(js4, [S => 999], (0, 1000.0), [β => 100.0, γ => 0.01])
jprob = JumpProblem(js4, dprob, Direct(), rng = rng)
m4 = getmean(jprob, Nsims)
@test abs(m4 - 2.0 / 0.01) * 0.01 / 2.0 < 0.01

# test second order rx runs
maj1 = MassActionJump(2.0, [0 => 1], [S => 1])
maj2 = MassActionJump(γ, [S => 2], [S => -1])
@named js4 = JumpSystem([maj1, maj2], t, [S], [β, γ])
js4 = complete(js4)
dprob = DiscreteProblem(js4, [S => 999], (0, 1000.0), [β => 100.0, γ => 0.01])
jprob = JumpProblem(js4, dprob, Direct(), rng = rng)
sol = solve(jprob, SSAStepper());

# issue #819
@testset "Combined system name collisions" begin
    sys1 = JumpSystem([maj1, maj2], t, [S], [β, γ], name = :sys1)
    sys2 = JumpSystem([maj1, maj2], t, [S], [β, γ], name = :sys1)
    @test_throws ArgumentError JumpSystem([sys1.γ ~ sys2.γ], t, [], [],
        systems = [sys1, sys2], name = :foo)
end

# test if param mapper is setup correctly for callbacks
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
jprob = JumpProblem(js5, dprob, Direct(), save_positions = (false, false), rng = rng)
@test all(jprob.massaction_jump.scaled_rates .== [1.0, 0.0])

pcondit(u, t, integrator) = t == 1000.0
function paffect!(integrator)
    integrator.ps[k1] = 0.0
    integrator.ps[k2] = 1.0
    reset_aggregated_jumps!(integrator)
end
sol = solve(jprob, SSAStepper(), tstops = [1000.0],
    callback = DiscreteCallback(pcondit, paffect!))
@test_skip sol.u[end][1] == 100 # TODO: Fix mass-action jumps in JumpProcesses

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
