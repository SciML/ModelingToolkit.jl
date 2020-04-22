using ModelingToolkit, DiffEqBase, DiffEqJump, Test, LinearAlgebra
MT = ModelingToolkit

# basic SIR model with tweaks
@parameters β γ t
@variables S I R
rate₁   = β*S*I
affect₁ = [S ~ S - 1, I ~ I + 1]
rate₂   = γ*I+t
affect₂ = [I ~ I - 1, R ~ R + 1]
j₁      = ConstantRateJump(rate₁,affect₁)
j₂      = VariableRateJump(rate₂,affect₂)
js      = JumpSystem([j₁,j₂], t, [S,I,R], [β,γ])
mtjump1 = MT.assemble_crj(js, j₁)
mtjump2 = MT.assemble_vrj(js, j₂)

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
mutable struct TestInt
    u
    p
    t
end
mtintegrator = TestInt(u,p,tf)
integrator   = TestInt(u,p,tf)
@test abs(mtjump1.rate(u,p,tf) - jump1.rate(u,p,tf)) < 10*eps()
@test abs(mtjump2.rate(u,p,tf) - jump2.rate(u,p,tf)) < 10*eps()
mtjump1.affect!(mtintegrator)
jump1.affect!(integrator)
@test norm(integrator.u - mtintegrator.u) < 10*eps()
mtintegrator.u .= u; integrator.u .= u
mtjump2.affect!(mtintegrator)
jump2.affect!(integrator)
@test norm(integrator.u - mtintegrator.u) < 10*eps()


# test can make and solve a jump problem
rate₂   = γ*I
affect₂ = [I ~ I - 1, R ~ R + 1]
j₃ = ConstantRateJump(rate₂,affect₂)
js2 = JumpSystem([j₁,j₃], t, [S,I,R], [β,γ])
u₀ = [999,1,0]; p = (0.1/1000,0.01); tspan = (0.,250.)
dprob = DiscreteProblem(u₀,tspan,p)
jprob = JumpProblem(js2, dprob, Direct())
sol = solve(jprob, SSAStepper())

using Plots
plot(sol)