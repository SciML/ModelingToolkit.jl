using ModelingToolkit, Unitful, OrdinaryDiffEq, DiffEqJump
using Test
MT = ModelingToolkit
@parameters τ [unit = u"ms"]
@variables t [unit = u"ms"] E(t) [unit = u"kJ"] P(t) [unit = u"MW"]
D = Differential(t)

@test MT.get_units(t) == 1u"ms"
@test MT.get_units(E) == 1u"kJ"
@test MT.get_units(τ) == 1u"ms"

@test MT.get_units(0.5) == 1.0
@test MT.get_units(t) == 1.0u"ms"
@test MT.get_units(P) == 1.0u"MW"
@test MT.get_units(τ) == 1.0u"ms"

@test MT.get_units(τ^-1) == 1/u"ms"
@test MT.get_units(D(E)) == 1.0u"MW"
@test MT.get_units(E/τ) == 1.0u"MW"
@test MT.get_units(2*P) == 1.0u"MW"
@test MT.get_units(t/τ) == 1.0
@test MT.get_units(P - E/τ)/1.0u"MW" == 1.0

@test MT.get_units(1.0^(t/τ)) == 1.0
@test MT.get_units(exp(t/τ)) == 1.0
@test MT.get_units(sin(t/τ)) == 1.0
@test MT.get_units(sin(1u"rad")) == 1.0
@test MT.get_units(t^2) == 1.0u"ms"^2

@test !MT.validate(E^1.5 ~ E^(t/τ))
@test MT.validate(E^(t/τ) ~ E^(t/τ))

eqs = [D(E) ~ P - E/τ
        0.0u"MW" ~ P]
@test MT.get_units(eqs[1].lhs) == 1.0u"MW"
@test MT.get_units(eqs[1].rhs) == 1.0u"MW"
@test MT.validate(eqs[1])
@test MT.validate(eqs[2])
@test MT.validate(eqs)
sys = ODESystem(eqs)
sys = ODESystem(eqs, t, [P, E], [τ])

@test !MT.validate(D(D(E)) ~ P)
@test !MT.validate(0 ~ P + E*τ)
@test_logs (:warn,) MT.validate(0 ~ P + E*τ)
@test_logs (:warn,) MT.validate(P + E*τ ~ 0)
@test_logs (:warn,) MT.validate(P ~ 0)

#Unit-free
@variables x y z u
@parameters σ ρ β
eqs = [0 ~ σ*(y - x)]
@test MT.validate(eqs)

#Array variables
@variables t x[1:3,1:3](t)
D = Differential(t)
eqs = D.(x) .~ x
ODESystem(eqs)

# Array ops
using Symbolics: unwrap, wrap
using LinearAlgebra
@variables t
sts = @variables x[1:3](t) y(t)
ps = @parameters p[1:3] = [1, 2, 3]
D = Differential(t)
eqs = [
       collect(D.(x) ~ x)
       D(y) ~ norm(x)*y
      ]
ODESystem(eqs, t, [sts...;], [ps...;])

#= Not supported yet b/c iterate doesn't work on unitful array
# Array ops with units
@variables t [unit =u"s"]
sts = @variables x[1:3](t) [unit = u"kg"] y(t) [unit = u"kg"]
ps = @parameters b [unit = u"s"^-1]
D = Differential(t)
eqs = [
       collect(D.(x) ~ b*x)
       D(y) ~ b*norm(x)
      ]
ODESystem(eqs, t, [sts...;], [ps...;])

#Array variables with units
@variables t [unit = u"s"] x[1:3,1:3](t) [unit = u"kg"] 
@parameters a [unit = u"s"^-1]
D = Differential(t)
eqs = D.(x) .~ a*x
ODESystem(eqs)
=#

#Difference equation with units
@parameters t [unit = u"s"] a [unit = u"s"^-1]
@variables x(t) [unit = u"kg"]
δ = Differential(t)
D = Difference(t; dt = 0.1u"s")
eqs = [
    δ(x) ~ a*x 
]
de = ODESystem(eqs, t, [x, y], [a])


@parameters t
@variables y[1:3](t)
@parameters k[1:3]
D = Differential(t)

eqs = [D(y[1]) ~ -k[1]*y[1] + k[3]*y[2]*y[3],
       D(y[2]) ~  k[1]*y[1] - k[3]*y[2]*y[3] - k[2]*y[2]^2,
       0 ~  y[1] + y[2] + y[3] - 1]

sys = ODESystem(eqs,t,y,k)

# Nonlinear system
@parameters a [unit = u"kg"^-1]
@variables x [unit = u"kg"]
eqs = [
    0 ~ a*x 
]
nls = NonlinearSystem(eqs, [x], [a])

# SDE test w/ noise vector
@parameters τ [unit = u"ms"] Q [unit = u"MW"]
@variables t [unit = u"ms"] E(t) [unit = u"kJ"] P(t) [unit = u"MW"]
D = Differential(t)
eqs = [D(E) ~ P - E/τ
      P ~ Q]

noiseeqs = [0.1u"MW",
            0.1u"MW"]
sys = SDESystem(eqs, noiseeqs, t, [P, E], [τ, Q])
# With noise matrix
noiseeqs = [0.1u"MW" 0.1u"MW"
            0.1u"MW" 0.1u"MW"]
sys = SDESystem(eqs,noiseeqs, t, [P, E], [τ, Q])

noiseeqs = [0.1u"MW" 0.1u"MW"
            0.1u"MW" 0.1u"s"]
@test !MT.validate(eqs,noiseeqs)

#Test non-trivial simplifications
@variables t [unit = u"s"] V(t) [unit = u"m"^3] L(t) [unit = u"m"]
@parameters v [unit = u"m/s"] r [unit =u"m"^3/u"s"]
D = Differential(t)
eqs = [D(L) ~ v,
       V ~ L^3]
sys = ODESystem(eqs)
sys_simple = structural_simplify(sys)

eqs = [D(V) ~ r,
       V ~ L^3]
sys = ODESystem(eqs)
sys_simple = structural_simplify(sys)

@variables V [unit = u"m"^3] L [unit = u"m"]
@parameters v [unit = u"m/s"] r [unit =u"m"^3/u"s"] t [unit = u"s"]
eqs = [V ~ r*t,
       V ~ L^3]
sys = NonlinearSystem(eqs, [V, L], [t, r])
sys_simple = structural_simplify(sys)

eqs = [L ~ v*t,
       V ~ L^3]
sys = NonlinearSystem(eqs, [V,L], [t,r])
sys_simple = structural_simplify(sys)

#Jump System
@parameters β [unit = 1/(u"mol"^2*u"s")] γ [unit = 1/(u"mol"*u"s")] t [unit = u"s"] 
@variables S(t) [unit = u"mol"] I(t) [unit = u"mol"] R(t) [unit = u"mol"]
rate₁   = β*S*I
affect₁ = [S ~ S - 1u"mol", I ~ I + 1u"mol"]
rate₂   = γ*I
affect₂ = [I ~ I - 1u"mol", R ~ R + 1u"mol"]
j₁      = ConstantRateJump(rate₁, affect₁)
j₂      = VariableRateJump(rate₂, affect₂)
js      = JumpSystem([j₁, j₂], t, [S, I, R], [β, γ])

affect_wrong = [S ~ S - 1u"mol", I ~ I + 1]
j_wrong      = ConstantRateJump(rate₁, affect_wrong)
@test_throws ArgumentError JumpSystem([j_wrong, j₂], t, [S, I, R], [β, γ])

rate_wrong   = γ^2*I
j_wrong     = ConstantRateJump(rate_wrong, affect₂)
@test_throws ArgumentError JumpSystem([j₁, j_wrong], t, [S, I, R], [β, γ])

# mass action jump tests for SIR model
maj1 = MassActionJump(2*β/2, [S => 1, I => 1], [S => -1, I => 1])
maj2 = MassActionJump(γ, [I => 1], [I => -1, R => 1])
js3  = JumpSystem([maj1, maj2], t, [S,I,R], [β,γ])

#Test unusual jump system
@parameters β γ t
@variables S(t) I(t) R(t)

maj1 = MassActionJump(2.0, [0 => 1], [S => 1])
maj2 = MassActionJump(γ, [S => 1], [S => -1])
js4  = JumpSystem([maj1, maj2], t, [S], [β, γ])