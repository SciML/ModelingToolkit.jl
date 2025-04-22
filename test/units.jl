using ModelingToolkit, OrdinaryDiffEq, JumpProcesses, Unitful
using Test
MT = ModelingToolkit
UMT = ModelingToolkit.UnitfulUnitCheck
@independent_variables t [unit = u"ms"]
@parameters τ [unit = u"ms"] γ
@variables E(t) [unit = u"kJ"] P(t) [unit = u"MW"]
D = Differential(t)

#This is how equivalent works:
@test UMT.equivalent(u"MW", u"kJ/ms")
@test !UMT.equivalent(u"m", u"cm")
@test UMT.equivalent(UMT.get_unit(P^γ), UMT.get_unit((E / τ)^γ))

# Basic access
@test UMT.get_unit(t) == u"ms"
@test UMT.get_unit(E) == u"kJ"
@test UMT.get_unit(τ) == u"ms"
@test UMT.get_unit(γ) == UMT.unitless
@test UMT.get_unit(0.5) == UMT.unitless
@test UMT.get_unit(UMT.SciMLBase.NullParameters()) == UMT.unitless

# Prohibited unit types
@parameters β [unit = u"°"] α [unit = u"°C"] γ [unit = 1u"s"]
@test_throws UMT.ValidationError UMT.get_unit(β)
@test_throws UMT.ValidationError UMT.get_unit(α)
@test_throws UMT.ValidationError UMT.get_unit(γ)

# Non-trivial equivalence & operators
@test UMT.get_unit(τ^-1) == u"ms^-1"
@test UMT.equivalent(UMT.get_unit(D(E)), u"MW")
@test UMT.equivalent(UMT.get_unit(E / τ), u"MW")
@test UMT.get_unit(2 * P) == u"MW"
@test UMT.get_unit(t / τ) == UMT.unitless
@test UMT.equivalent(UMT.get_unit(P - E / τ), u"MW")
@test UMT.equivalent(UMT.get_unit(D(D(E))), u"MW/ms")
@test UMT.get_unit(ifelse(t > t, P, E / τ)) == u"MW"
@test UMT.get_unit(1.0^(t / τ)) == UMT.unitless
@test UMT.get_unit(exp(t / τ)) == UMT.unitless
@test UMT.get_unit(sin(t / τ)) == UMT.unitless
@test UMT.get_unit(sin(1 * u"rad")) == UMT.unitless
@test UMT.get_unit(t^2) == u"ms^2"

eqs = [D(E) ~ P - E / τ
       0 ~ P]
@test UMT.validate(eqs)
@named sys = ODESystem(eqs, t)

@test !UMT.validate(D(D(E)) ~ P)
@test !UMT.validate(0 ~ P + E * τ)

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
        x(t)=1.0, [unit = NoUnits])
    ODESystem(Equation[], t, sts, []; name = name)
end
@named p1 = Pin()
@named p2 = Pin()
@named op = OtherPin()
@named lp = LongPin()
good_eqs = [connect(p1, p2)]
bad_eqs = [connect(p1, p2, op)]
bad_length_eqs = [connect(op, lp)]
@test UMT.validate(good_eqs)
@test !UMT.validate(bad_eqs)
@test !UMT.validate(bad_length_eqs)
@named sys = ODESystem(good_eqs, t, [], [])
@test_throws MT.ValidationError ODESystem(bad_eqs, t, [], []; name = :sys)

# Array variables
@independent_variables t [unit = u"s"]
@parameters v[1:3]=[1, 2, 3] [unit = u"m/s"]
@variables x(t)[1:3] [unit = u"m"]
D = Differential(t)
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
@independent_variables t [unit = u"ms"]
@parameters τ [unit = u"ms"] Q [unit = u"MW"]
@variables E(t) [unit = u"kJ"] P(t) [unit = u"MW"]
D = Differential(t)
eqs = [D(E) ~ P - E / τ
       P ~ Q]

noiseeqs = [0.1u"MW",
    0.1u"MW"]
@named sys = SDESystem(eqs, noiseeqs, t, [P, E], [τ, Q])

# With noise matrix
noiseeqs = [0.1u"MW" 0.1u"MW"
            0.1u"MW" 0.1u"MW"]
@named sys = SDESystem(eqs, noiseeqs, t, [P, E], [τ, Q])

# Invalid noise matrix
noiseeqs = [0.1u"MW" 0.1u"MW"
            0.1u"MW" 0.1u"s"]
@test !UMT.validate(eqs, noiseeqs)

# Non-trivial simplifications
@independent_variables t [unit = u"s"]
@parameters v [unit = u"m/s"] r [unit = u"m"^3 / u"s"]
@variables V(t) [unit = u"m"^3] L(t) [unit = u"m"]
D = Differential(t)
eqs = [D(L) ~ v,
    V ~ L^3]
@named sys = ODESystem(eqs, t)
sys_simple = structural_simplify(sys)

eqs = [D(V) ~ r,
    V ~ L^3]
@named sys = ODESystem(eqs, t)
sys_simple = structural_simplify(sys)

@variables V [unit = u"m"^3] L [unit = u"m"]
@parameters v [unit = u"m/s"] r [unit = u"m"^3 / u"s"] t [unit = u"s"]
eqs = [V ~ r * t,
    V ~ L^3]
@named sys = NonlinearSystem(eqs, [V, L], [t, r])
sys_simple = structural_simplify(sys)

eqs = [L ~ v * t,
    V ~ L^3]
@named sys = NonlinearSystem(eqs, [V, L], [t, r])
sys_simple = structural_simplify(sys)

#Jump System
@parameters β [unit = u"(mol^2*s)^-1"] γ [unit = u"(mol*s)^-1"] t [unit = u"s"] jumpmol [
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
@parameters β γ t
@variables S(t) I(t) R(t)

maj1 = MassActionJump(2.0, [0 => 1], [S => 1])
maj2 = MassActionJump(γ, [S => 1], [S => -1])
@named js4 = JumpSystem([maj1, maj2], t, [S], [β, γ])

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

@variables x(t)
@test ModelingToolkit.get_unit(sin(x)) == ModelingToolkit.unitless

@mtkmodel ExpressionParametersTest begin
    @parameters begin
        v = 1.0, [unit = u"m/s"]
        τ = 1.0, [unit = u"s"]
    end
    @components begin
        pt = ParamTest(; a = v * τ)
    end
end

@named sys = ExpressionParametersTest(; v = 2.0u"m/s", τ = 3.0u"s")
sys = complete(sys)
# TODO: Is there a way to evaluate this expression and compare to 6.0?
@test isequal(ModelingToolkit.getdefault(sys.pt.a), sys.v * sys.τ)
@test ModelingToolkit.getdefault(sys.v) ≈ 2.0
@test ModelingToolkit.getdefault(sys.τ) ≈ 3.0
