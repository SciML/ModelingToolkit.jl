using ModelingToolkit, DiffEqJump, Test, LinearAlgebra
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
mtjump2 = MT.assemble_crj(js, j₂)

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
t = 1.0
mutable struct TestInt
    u
    p
    t
end
mtintegrator = TestInt(u,p,t)
integrator   = TestInt(u,p,t)
@test abs(mtjump1.rate(u,p,t) - jump1.rate(u,p,t)) < 10*eps()
@test abs(mtjump2.rate(u,p,t) - jump2.rate(u,p,t)) < 10*eps()
mtjump1.affect!(mtintegrator)
jump1.affect!(integrator)
@test norm(integrator.u - mtintegrator.u) < 10*eps()
mtintegrator.u .= u; integrator.u .= u
mtjump2.affect!(mtintegrator)
jump2.affect!(integrator)
@test norm(integrator.u - mtintegrator.u) < 10*eps()


