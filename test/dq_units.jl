using ModelingToolkit, OrdinaryDiffEq, JumpProcesses, DynamicQuantities
using Test
MT = ModelingToolkit
using ModelingToolkit: t, D
@parameters τ [unit = u"s"] γ
@variables E(t) [unit = u"J"] P(t) [unit = u"W"]

# Basic access
@test MT.get_unit(t) == u"s"
@test MT.get_unit(E) == u"J"
@test MT.get_unit(τ) == u"s"
@test MT.get_unit(γ) == MT.unitless
@test MT.get_unit(0.5) == MT.unitless
@test MT.get_unit(MT.SciMLBase.NullParameters()) == MT.unitless

# Prohibited unit types
@parameters γ [unit = 1u"ms"]
@test_throws MT.ValidationError MT.get_unit(γ)

eqs = [D(E) ~ P - E / τ
       0 ~ P]
@test MT.validate(eqs)
@named sys = ODESystem(eqs, t)

@test !MT.validate(D(D(E)) ~ P)
@test !MT.validate(0 ~ P + E * τ)

# Disabling unit validation/checks selectively
@test_throws MT.ArgumentError ODESystem(eqs, t, [E, P, t], [τ], name = :sys)
ODESystem(eqs, t, [E, P, t], [τ], name = :sys, checks = MT.CheckUnits)
eqs = [D(E) ~ P - E / τ
       0 ~ P + E * τ]
@test_throws MT.ValidationError ODESystem(eqs, t, name = :sys, checks = MT.CheckAll)
@test_throws MT.ValidationError ODESystem(eqs, t, name = :sys, checks = true)
ODESystem(eqs, t, name = :sys, checks = MT.CheckNone)
ODESystem(eqs, t, name = :sys, checks = false)
@test_throws MT.ValidationError ODESystem(eqs, t, name = :sys,
    checks = MT.CheckComponents | MT.CheckUnits)
@named sys = ODESystem(eqs, t, checks = MT.CheckComponents)
@test_throws MT.ValidationError ODESystem(eqs, t, [E, P, t], [τ], name = :sys,
    checks = MT.CheckUnits)

# connection validation
@connector function Pin(; name)
    sts = @variables(v(t)=1.0, [unit = u"V"],
        i(t)=1.0, [unit = u"A", connect = Flow])
    ODESystem(Equation[], t, sts, []; name = name)
end
@connector function OtherPin(; name)
    sts = @variables(v(t)=1.0, [unit = u"mV"],
        i(t)=1.0, [unit = u"mA", connect = Flow])
    ODESystem(Equation[], t, sts, []; name = name)
end
@connector function LongPin(; name)
    sts = @variables(v(t)=1.0, [unit = u"V"],
        i(t)=1.0, [unit = u"A", connect = Flow],
        x(t)=1.0)
    ODESystem(Equation[], t, sts, []; name = name)
end
@named p1 = Pin()
@named p2 = Pin()
@test_throws MT.ValidationError @named op = OtherPin()
@named lp = LongPin()
good_eqs = [connect(p1, p2)]
@test MT.validate(good_eqs)
@named sys = ODESystem(good_eqs, t, [], [])

# Array variables
@variables x(t)[1:3] [unit = u"m"]
@parameters v[1:3]=[1, 2, 3] [unit = u"m/s"]
eqs = D.(x) .~ v
ODESystem(eqs, t, name = :sys)

# Nonlinear system
@parameters a [unit = u"kg"^-1]
@variables x [unit = u"kg"]
eqs = [
    0 ~ a * x
]
@named nls = NonlinearSystem(eqs, [x], [a])

# SDE test w/ noise vector
@parameters τ [unit = u"s"] Q [unit = u"W"]
@variables E(t) [unit = u"J"] P(t) [unit = u"W"]
eqs = [D(E) ~ P - E / τ
       P ~ Q]

noiseeqs = [0.1u"W",
    0.1u"W"]
@named sys = SDESystem(eqs, noiseeqs, t, [P, E], [τ, Q])

# With noise matrix
noiseeqs = [0.1u"W" 0.1u"W"
            0.1u"W" 0.1u"W"]
@named sys = SDESystem(eqs, noiseeqs, t, [P, E], [τ, Q])

# Invalid noise matrix
noiseeqs = [0.1u"W" 0.1u"W"
            0.1u"W" 0.1u"s"]
@test !MT.validate(eqs, noiseeqs)

# Non-trivial simplifications
@variables V(t) [unit = u"m"^3] L(t) [unit = u"m"]
@parameters v [unit = u"m/s"] r [unit = u"m"^3 / u"s"]
eqs = [D(L) ~ v,
    V ~ L^3]
@named sys = ODESystem(eqs, t)
sys_simple = structural_simplify(sys)

eqs = [D(V) ~ r,
    V ~ L^3]
@named sys = ODESystem(eqs, t)
sys_simple = structural_simplify(sys)

@variables V [unit = u"m"^3] L [unit = u"m"]
@parameters v [unit = u"m/s"] r [unit = u"m"^3 / u"s"]
eqs = [V ~ r * t,
    V ~ L^3]
@named sys = NonlinearSystem(eqs, [V, L], [t, r])
sys_simple = structural_simplify(sys)

eqs = [L ~ v * t,
    V ~ L^3]
@named sys = NonlinearSystem(eqs, [V, L], [t, r])
sys_simple = structural_simplify(sys)

#Jump System
@parameters β [unit = u"(mol^2*s)^-1"] γ [unit = u"(mol*s)^-1"] jumpmol [
    unit = u"mol"
]
@variables S(t) [unit = u"mol"] I(t) [unit = u"mol"] R(t) [unit = u"mol"]
rate₁ = β * S * I
affect₁ = [S ~ S - 1 * jumpmol, I ~ I + 1 * jumpmol]
rate₂ = γ * I
affect₂ = [I ~ I - 1 * jumpmol, R ~ R + 1 * jumpmol]
j₁ = ConstantRateJump(rate₁, affect₁)
j₂ = VariableRateJump(rate₂, affect₂)
js = JumpSystem([j₁, j₂], t, [S, I, R], [β, γ], name = :sys)

affect_wrong = [S ~ S - jumpmol, I ~ I + 1]
j_wrong = ConstantRateJump(rate₁, affect_wrong)
@test_throws MT.ValidationError JumpSystem([j_wrong, j₂], t, [S, I, R], [β, γ], name = :sys)

rate_wrong = γ^2 * I
j_wrong = ConstantRateJump(rate_wrong, affect₂)
@test_throws MT.ValidationError JumpSystem([j₁, j_wrong], t, [S, I, R], [β, γ], name = :sys)

# mass action jump tests for SIR model
maj1 = MassActionJump(2 * β / 2, [S => 1, I => 1], [S => -1, I => 1])
maj2 = MassActionJump(γ, [I => 1], [I => -1, R => 1])
@named js3 = JumpSystem([maj1, maj2], t, [S, I, R], [β, γ])

#Test unusual jump system
@parameters β γ
@variables S(t) I(t) R(t)

maj1 = MassActionJump(2.0, [0 => 1], [S => 1])
maj2 = MassActionJump(γ, [S => 1], [S => -1])
@named js4 = JumpSystem([maj1, maj2], ModelingToolkit.t_nounits, [S], [β, γ])

@mtkmodel ParamTest begin
    @parameters begin
        a, [unit = u"m"]
    end
    @variables begin
        b(t), [unit = u"kg"]
    end
end

@named sys = ParamTest()

@named sys = ParamTest(a = 3.0u"cm")
@test ModelingToolkit.getdefault(sys.a) ≈ 0.03

@test_throws ErrorException ParamTest(; name = :t, a = 1.0)
@test_throws ErrorException ParamTest(; name = :t, a = 1.0u"s")

@mtkmodel ArrayParamTest begin
    @parameters begin
        a[1:2], [unit = u"m"]
    end
end

@named sys = ArrayParamTest()

@named sys = ArrayParamTest(a = [1.0, 3.0]u"cm")
@test ModelingToolkit.getdefault(sys.a) ≈ [0.01, 0.03]
