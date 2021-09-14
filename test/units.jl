using ModelingToolkit, Unitful, OrdinaryDiffEq, DiffEqJump, IfElse
using Test
MT = ModelingToolkit
@parameters τ [unit = u"ms"] γ
@variables t [unit = u"ms"] E(t) [unit = u"kJ"] P(t) [unit = u"MW"]
D = Differential(t)

#This is how equivalent works:
@test MT.equivalent(u"MW" ,u"kJ/ms")
@test !MT.equivalent(u"m", u"cm")
@test MT.equivalent(MT.get_unit(P^γ), MT.get_unit((E/τ)^γ))

# Basic access
@test MT.get_unit(t) == u"ms"
@test MT.get_unit(E) == u"kJ"
@test MT.get_unit(τ) == u"ms"
@test MT.get_unit(γ) == MT.unitless
@test MT.get_unit(0.5) == MT.unitless
@test MT.get_unit(MT.SciMLBase.NullParameters) == MT.unitless

# Prohibited unit types
@parameters β [unit = u"°"] α [unit = u"°C"] γ [unit = 1u"s"]
@test_throws MT.ValidationError MT.get_unit(β)
@test_throws MT.ValidationError MT.get_unit(α)
@test_throws MT.ValidationError MT.get_unit(γ)

# Non-trivial equivalence & operators
@test MT.get_unit(τ^-1) == u"ms^-1"
@test MT.equivalent(MT.get_unit(D(E)),u"MW")
@test MT.equivalent(MT.get_unit(E/τ), u"MW")
@test MT.get_unit(2*P) == u"MW"
@test MT.get_unit(t/τ) == MT.unitless
@test MT.equivalent(MT.get_unit(P - E/τ),u"MW")
@test MT.equivalent(MT.get_unit(D(D(E))),u"MW/ms")
@test MT.get_unit(IfElse.ifelse(t>t,P,E/τ)) == u"MW"
@test MT.get_unit(1.0^(t/τ)) == MT.unitless
@test MT.get_unit(exp(t/τ)) == MT.unitless
@test MT.get_unit(sin(t/τ)) == MT.unitless
@test MT.get_unit(sin(1u"rad")) == MT.unitless
@test MT.get_unit(t^2) == u"ms^2"

eqs = [D(E) ~ P - E/τ
        0 ~ P]
@test MT.validate(eqs)
@named sys = ODESystem(eqs)

@test !MT.validate(D(D(E)) ~ P)
@test !MT.validate(0 ~ P + E*τ)

# Array variables
@variables t [unit = u"s"] x[1:3](t) [unit = u"m"]
@parameters v[1:3] = [1,2,3] [unit = u"m/s"]
D = Differential(t)
eqs = D.(x) .~ v
ODESystem(eqs,name=:sys)

# Difference equation
@parameters t [unit = u"s"] a [unit = u"s"^-1]
@variables x(t) [unit = u"kg"]
δ = Differential(t)
D = Difference(t; dt = 0.1u"s")
eqs = [
    δ(x) ~ a*x 
]
de = ODESystem(eqs, t, [x], [a],name=:sys)

# Nonlinear system
@parameters a [unit = u"kg"^-1]
@variables x [unit = u"kg"]
eqs = [
    0 ~ a*x 
]
@named nls = NonlinearSystem(eqs, [x], [a])

# SDE test w/ noise vector
@parameters τ [unit = u"ms"] Q [unit = u"MW"]
@variables t [unit = u"ms"] E(t) [unit = u"kJ"] P(t) [unit = u"MW"]
D = Differential(t)
eqs = [D(E) ~ P - E/τ
      P ~ Q]

noiseeqs = [0.1u"MW",
            0.1u"MW"]
@named sys = SDESystem(eqs, noiseeqs, t, [P, E], [τ, Q])

# With noise matrix
noiseeqs = [0.1u"MW" 0.1u"MW"
            0.1u"MW" 0.1u"MW"]
@named sys = SDESystem(eqs,noiseeqs, t, [P, E], [τ, Q])

# Invalid noise matrix 
noiseeqs = [0.1u"MW" 0.1u"MW"
            0.1u"MW" 0.1u"s"]
@test !MT.validate(eqs,noiseeqs)

# Non-trivial simplifications
@variables t [unit = u"s"] V(t) [unit = u"m"^3] L(t) [unit = u"m"]
@parameters v [unit = u"m/s"] r [unit =u"m"^3/u"s"]
D = Differential(t)
eqs = [D(L) ~ v,
       V ~ L^3]
@named sys = ODESystem(eqs)
sys_simple = structural_simplify(sys)

eqs = [D(V) ~ r,
       V ~ L^3]
@named sys = ODESystem(eqs)
sys_simple = structural_simplify(sys)

@variables V [unit = u"m"^3] L [unit = u"m"]
@parameters v [unit = u"m/s"] r [unit = u"m"^3/u"s"] t [unit = u"s"]
eqs = [V ~ r*t,
       V ~ L^3]
@named sys = NonlinearSystem(eqs, [V, L], [t, r])
sys_simple = structural_simplify(sys)

eqs = [L ~ v*t,
       V ~ L^3]
@named sys = NonlinearSystem(eqs, [V,L], [t,r])
sys_simple = structural_simplify(sys)

#Jump System
@parameters β [unit = u"(mol^2*s)^-1"] γ [unit = u"(mol*s)^-1"] t [unit = u"s"] jumpmol [unit = u"mol"]
@variables S(t) [unit = u"mol"] I(t) [unit = u"mol"] R(t) [unit = u"mol"]
rate₁   = β*S*I
affect₁ = [S ~ S - 1*jumpmol, I ~ I + 1*jumpmol]
rate₂   = γ*I
affect₂ = [I ~ I - 1*jumpmol, R ~ R + 1*jumpmol]
j₁      = ConstantRateJump(rate₁, affect₁)
j₂      = VariableRateJump(rate₂, affect₂)
js      = JumpSystem([j₁, j₂], t, [S, I, R], [β, γ],name=:sys)

affect_wrong = [S ~ S - jumpmol, I ~ I + 1]
j_wrong      = ConstantRateJump(rate₁, affect_wrong)
@test_throws MT.ValidationError JumpSystem([j_wrong, j₂], t, [S, I, R], [β, γ],name=:sys)

rate_wrong   = γ^2*I
j_wrong     = ConstantRateJump(rate_wrong, affect₂)
@test_throws MT.ValidationError JumpSystem([j₁, j_wrong], t, [S, I, R], [β, γ],name=:sys)

# mass action jump tests for SIR model
maj1 = MassActionJump(2*β/2, [S => 1, I => 1], [S => -1, I => 1])
maj2 = MassActionJump(γ, [I => 1], [I => -1, R => 1])
@named js3  = JumpSystem([maj1, maj2], t, [S,I,R], [β,γ])

#Test unusual jump system
@parameters β γ t
@variables S(t) I(t) R(t)

maj1 = MassActionJump(2.0, [0 => 1], [S => 1])
maj2 = MassActionJump(γ, [S => 1], [S => -1])
@named js4  = JumpSystem([maj1, maj2], t, [S], [β, γ])

