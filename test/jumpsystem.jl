using ModelingToolkit, DiffEqBase, DiffEqJump, Test, LinearAlgebra
MT = ModelingToolkit

# basic MT SIR model with tweaks
@parameters β γ t
@variables S I R
rate₁   = β*S*I
affect₁ = [S ~ S - 1, I ~ I + 1]
rate₂   = γ*I+t
affect₂ = [I ~ I - 1, R ~ R + 1]
j₁      = ConstantRateJump(rate₁,affect₁)
j₂      = VariableRateJump(rate₂,affect₂)
js      = JumpSystem([j₁,j₂], t, [S,I,R], [β,γ])
statetoid = Dict(convert(Variable,state) => i for (i,state) in enumerate(states(js)))
mtjump1 = MT.assemble_crj(js, j₁, statetoid)
mtjump2 = MT.assemble_vrj(js, j₂, statetoid)

# doc version
rate1(u,p,t) = (0.1/1000.0)*u[1]*u[2]
function affect1!(integrator)
  integrator.u[1] -= 1
  integrator.u[2] += 1
end
jump1 = ConstantRateJump(rate1,affect1!)
rate2(u,p,t) = 0.01u[2]+t
function affect2!(integrator)
  integrator.u[2] -= 1
  integrator.u[3] += 1
end
jump2 = VariableRateJump(rate2,affect2!)

# test crjs
u = [100, 9, 5]
p = (0.1/1000,0.01)
tf = 1.0
mutable struct TestInt{U,V,T}
    u::U
    p::V
    t::T
end
mtintegrator = TestInt(u,p,tf)
integrator   = TestInt(u,p,tf)
@test abs(mtjump1.rate(u,p,tf) - jump1.rate(u,p,tf)) < 10*eps()
@test abs(mtjump2.rate(u,p,tf) - jump2.rate(u,p,tf)) < 10*eps()
mtjump1.affect!(mtintegrator)
jump1.affect!(integrator)
@test all(integrator.u .== mtintegrator.u) 
mtintegrator.u .= u; integrator.u .= u
mtjump2.affect!(mtintegrator)
jump2.affect!(integrator)
@test all(integrator.u .== mtintegrator.u)

# test MT can make and solve a jump problem
rate₃   = γ*I
affect₃ = [I ~ I - 1, R ~ R + 1]
j₃ = ConstantRateJump(rate₃,affect₃)
js2 = JumpSystem([j₁,j₃], t, [S,I,R], [β,γ])
u₀ = [999,1,0]; p = (0.1/1000,0.01); tspan = (0.,250.)
u₀map = [S => 999, I => 1, R => 0]
parammap = [β => .1/1000, γ => .01]
dprob = DiscreteProblem(js2, u₀map, tspan, parammap)
jprob = JumpProblem(js2, dprob, Direct(), save_positions=(false,false))
Nsims = 10000
function getmean(jprob,Nsims)
  m = 0.0
  for i = 1:Nsims
    sol = solve(jprob, SSAStepper())
    m += sol[end,end]
  end
  m/Nsims
end
m = getmean(jprob,Nsims)

#test the MT JumpProblem rates/affects are correct
rate2(u,p,t) = 0.01u[2]
jump2 = ConstantRateJump(rate2,affect2!)
mtjumps = jprob.discrete_jump_aggregation
@test abs(mtjumps.rates[1](u,p,tf) - jump1.rate(u,p,tf)) < 10*eps()
@test abs(mtjumps.rates[2](u,p,tf) - jump2.rate(u,p,tf)) < 10*eps()
mtjumps.affects![1](mtintegrator)
jump1.affect!(integrator)
@test all(integrator.u .== mtintegrator.u) 
mtintegrator.u .= u; integrator.u .= u
mtjumps.affects![2](mtintegrator)
jump2.affect!(integrator)
@test all(integrator.u .== mtintegrator.u)

# direct vers
p = (0.1/1000,0.01)
prob = DiscreteProblem([999,1,0],(0.0,250.0),p)
r1(u,p,t) = (0.1/1000.0)*u[1]*u[2]
function a1!(integrator)
  integrator.u[1] -= 1
  integrator.u[2] += 1
end
j1 = ConstantRateJump(r1,a1!)
r2(u,p,t) = 0.01u[2]
function a2!(integrator)
  integrator.u[2] -= 1
  integrator.u[3] += 1
end
j2 = ConstantRateJump(r2,a2!)
jset = JumpSet((),(j1,j2),nothing,nothing)
jprob = JumpProblem(prob,Direct(),jset, save_positions=(false,false))
m2 = getmean(jprob,Nsims)

# test JumpSystem solution agrees with direct version
@test abs(m-m2)/m < .01

# mass action jump tests for SIR model
maj1 = MassActionJump(2*β/2, [S => 1, I => 1], [S => -1, I => 1])
maj2 = MassActionJump(γ, [I => 1], [I => -1, R => 1])
js3   = JumpSystem([maj1,maj2], t, [S,I,R], [β,γ])
statetoid = Dict(convert(Variable,state) => i for (i,state) in enumerate(states(js)))
ptoid     = Dict(convert(Variable,par) => i for (i,par) in enumerate(parameters(js)))
dprob = DiscreteProblem(js3, u₀map, tspan, parammap)
jprob = JumpProblem(js3, dprob, Direct())
m3 = getmean(jprob,Nsims)
@test abs(m-m3)/m < .01

# mass action jump tests for other reaction types (zero order, decay)
maj1 = MassActionJump(2.0, [0 => 1], [S => 1])
maj2 = MassActionJump(γ, [S => 1], [S => -1])
js4   = JumpSystem([maj1,maj2], t, [S], [β,γ])
statetoid = Dict(convert(Variable,state) => i for (i,state) in enumerate(states(js)))
ptoid     = Dict(convert(Variable,par) => i for (i,par) in enumerate(parameters(js)))
dprob = DiscreteProblem(js4, [S => 999], (0,1000.), [β => 100.,γ => .01])
jprob = JumpProblem(js4, dprob, Direct())
m4 = getmean(jprob,Nsims)
@test abs(m4 - 2.0/.01)*.01/2.0 < .01

# test second order rx runs 
maj1 = MassActionJump(2.0, [0 => 1], [S => 1])
maj2 = MassActionJump(γ, [S => 2], [S => -1])
js4   = JumpSystem([maj1,maj2], t, [S], [β,γ])
statetoid = Dict(convert(Variable,state) => i for (i,state) in enumerate(states(js)))
ptoid     = Dict(convert(Variable,par) => i for (i,par) in enumerate(parameters(js)))
dprob = DiscreteProblem(js4, [S => 999], (0,1000.), [β => 100.,γ => .01])
jprob = JumpProblem(js4, dprob, Direct())
sol = solve(jprob, SSAStepper())
