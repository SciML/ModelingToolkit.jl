using ModelingToolkit, DiffEqBase, JumpProcesses, Test, LinearAlgebra
MT = ModelingToolkit

# basic MT SIR model with tweaks
@parameters β γ t
@variables S(t) I(t) R(t)
rate₁ = β * S * I
affect₁ = [S ~ S - 1, I ~ I + 1]
rate₂ = γ * I + t
affect₂ = [I ~ I - 1, R ~ R + 1]
j₁ = ConstantRateJump(rate₁, affect₁)
j₂ = VariableRateJump(rate₂, affect₂)
@named js = JumpSystem([j₁, j₂], t, [S, I, R], [β, γ])
statetoid = Dict(MT.value(state) => i for (i, state) in enumerate(states(js)))
mtjump1 = MT.assemble_crj(js, j₁, statetoid)
mtjump2 = MT.assemble_vrj(js, j₂, statetoid)

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
rate₃ = γ * I
affect₃ = [I ~ I - 1, R ~ R + 1]
j₃ = ConstantRateJump(rate₃, affect₃)
@named js2 = JumpSystem([j₁, j₃], t, [S, I, R], [β, γ])
u₀ = [999, 1, 0];
p = (0.1 / 1000, 0.01);
tspan = (0.0, 250.0);
u₀map = [S => 999, I => 1, R => 0]
parammap = [β => 0.1 / 1000, γ => 0.01]
dprob = DiscreteProblem(js2, u₀map, tspan, parammap)
jprob = JumpProblem(js2, dprob, Direct(), save_positions = (false, false))
Nsims = 30000
function getmean(jprob, Nsims)
    m = 0.0
    for i in 1:Nsims
        sol = solve(jprob, SSAStepper())
        m += sol[end, end]
    end
    m / Nsims
end
m = getmean(jprob, Nsims)

@variables S2(t)
obs = [S2 ~ 2 * S]
@named js2b = JumpSystem([j₁, j₃], t, [S, I, R], [β, γ], observed = obs)
dprob = DiscreteProblem(js2b, u₀map, tspan, parammap)
jprob = JumpProblem(js2b, dprob, Direct(), save_positions = (false, false))
sol = solve(jprob, SSAStepper(), saveat = tspan[2] / 10)
@test all(2 .* sol[S] .== sol[S2])

# test save_positions is working
jprob = JumpProblem(js2, dprob, Direct(), save_positions = (false, false))
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
jprob = JumpProblem(prob, Direct(), jset, save_positions = (false, false))
m2 = getmean(jprob, Nsims)

# test JumpSystem solution agrees with direct version
@test abs(m - m2) / m < 0.01

# mass action jump tests for SIR model
maj1 = MassActionJump(2 * β / 2, [S => 1, I => 1], [S => -1, I => 1])
maj2 = MassActionJump(γ, [I => 1], [I => -1, R => 1])
@named js3 = JumpSystem([maj1, maj2], t, [S, I, R], [β, γ])
dprob = DiscreteProblem(js3, u₀map, tspan, parammap)
jprob = JumpProblem(js3, dprob, Direct())
m3 = getmean(jprob, Nsims)
@test abs(m - m3) / m < 0.01

# maj jump test with various dep graphs
@named js3b = JumpSystem([maj1, maj2], t, [S, I, R], [β, γ])
jprobb = JumpProblem(js3b, dprob, NRM())
m4 = getmean(jprobb, Nsims)
@test abs(m - m4) / m < 0.01
jprobc = JumpProblem(js3b, dprob, RSSA())
m4 = getmean(jprobc, Nsims)
@test abs(m - m4) / m < 0.01

# mass action jump tests for other reaction types (zero order, decay)
maj1 = MassActionJump(2.0, [0 => 1], [S => 1])
maj2 = MassActionJump(γ, [S => 1], [S => -1])
@named js4 = JumpSystem([maj1, maj2], t, [S], [β, γ])
dprob = DiscreteProblem(js4, [S => 999], (0, 1000.0), [β => 100.0, γ => 0.01])
jprob = JumpProblem(js4, dprob, Direct())
m4 = getmean(jprob, Nsims)
@test abs(m4 - 2.0 / 0.01) * 0.01 / 2.0 < 0.01

# test second order rx runs
maj1 = MassActionJump(2.0, [0 => 1], [S => 1])
maj2 = MassActionJump(γ, [S => 2], [S => -1])
@named js4 = JumpSystem([maj1, maj2], t, [S], [β, γ])
dprob = DiscreteProblem(js4, [S => 999], (0, 1000.0), [β => 100.0, γ => 0.01])
jprob = JumpProblem(js4, dprob, Direct())
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
p = [k1 => 2.0, k2 => 0.0, k3 => 0.5]
u₀ = [A => 100, B => 0]
tspan = (0.0, 2000.0)
dprob = DiscreteProblem(js5, u₀, tspan, p)
jprob = JumpProblem(js5, dprob, Direct(), save_positions = (false, false))
@test all(jprob.massaction_jump.scaled_rates .== [1.0, 0.0])

pcondit(u, t, integrator) = t == 1000.0
function paffect!(integrator)
    integrator.p[1] = 0.0
    integrator.p[2] = 1.0
    reset_aggregated_jumps!(integrator)
end
sol = solve(jprob, SSAStepper(), tstops = [1000.0],
            callback = DiscreteCallback(pcondit, paffect!))
@test sol[1, end] == 100

# observed variable handling
@variables OBS(t)
@named js5 = JumpSystem([maj1, maj2], t, [S], [β, γ]; observed = [OBS ~ 2 * S])
OBS2 = OBS
@test isequal(OBS2, @nonamespace js5.OBS)
@unpack OBS = js5
@test isequal(OBS2, OBS)
