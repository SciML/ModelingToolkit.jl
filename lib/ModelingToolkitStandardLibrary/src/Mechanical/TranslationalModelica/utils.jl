@connector Flange begin
    s(t)
    f(t), [connect = Flow]
end
Base.@doc """
    Flange(;name)

1-dim. translational flange.

# States:
- `s`: [m] Absolute position of flange
- `f`: [N] Cut force into the flange
""" Flange

@connector Support begin
    s(t)
    f(t), [connect = Flow]
end
Base.@doc """
    Support(;name)

Support/housing 1-dim. translational flange.

# States:
- `s`: [m] Absolute position of the support/housing
- `f`: [N] Cut force into the flange
""" Support

@mtkmodel PartialTwoFlanges begin
    @components begin
        flange_a = Flange() # (left) driving flange (flange axis directed into cut plane, e. g. from left to right)
        flange_b = Flange() # (right) driven flange (flange axis directed out of cut plane)
    end
end

"""
    PartialCompliant(; name, s_rel = 0.0, f = 0.0)

Partial model for the compliant connection of two translational 1-dim. flanges.

# States:

  - `s_rel`: [m] Relative distance (= flange_b.s - flange_a.s). It accepts an initial value, which defaults to 0.0.
  - `f`: [N] Force between flanges (= flange_b.f). It accepts an initial value, which defaults to 0.0.
"""
@mtkmodel PartialCompliant begin
    @extend (flange_a, flange_b) = pt = PartialTwoFlanges()
    @variables begin
        s_rel(t), [description = "Relative distance between flanges", guess = 0.0]
        f(t), [description = "Force between flanges", guess = 0.0]
    end

    @equations begin
        s_rel ~ flange_b.s - flange_a.s
        flange_b.f ~ +f
        flange_a.f ~ -f
    end
end

"""
    PartialCompliantWithRelativeStates(;name, s_rel = 0.0, v_rel = 0.0, f = 0.0)

Partial model for the compliant connection of two translational 1-dim. flanges.

    # States:

  - `s_rel`: [m] Relative distance (= flange_b.phi - flange_a.phi). It accepts an initial value, which defaults to 0.0.
  - `v_rel`: [m/s] Relative linear velocity (= der(s_rel)). It accepts an initial value, which defaults to 0.0.
  - `f`: [N] Force between flanges (= flange_b.f). It accepts an initial value, which defaults to 0.0.
"""
@mtkmodel PartialCompliantWithRelativeStates begin
    @extend flange_a, flange_b = pt = PartialTwoFlanges()
    @variables begin
        s_rel(t), [description = "Relative distance between flanges"]
        v_rel(t), [description = "Relative linear velocity))"]
        f(t), [description = "Forces between flanges"]
    end

    @equations begin
        s_rel ~ flange_b.s - flange_a.s
        v_rel ~ D(s_rel)
        flange_b.f ~ f
        flange_a.f ~ -f
    end
end

"""
    PartialElementaryOneFlangeAndSupport2(; name, use_support = false)

Partial model for a component with one translational 1-dim. shaft flange and a support used for textual modeling, i.e., for elementary models

# Parameters:

  - `use_support`: If support flange enabled, otherwise implicitly grounded

# States:

  - `s_support`: [m] Absolute position of support flange"
"""
function PartialElementaryOneFlangeAndSupport2(; name, use_support = false)
    @named flange = Flange()
    @variables s_support(t) [description = "Absolute position of support flange"]
    @variables s(t) [
        description = "Distance between flange and support (= flange.s - support.s)"
    ]
    eqs = [s ~ flange.s - s_support]
    if use_support
        @named support = Support()
        push!(eqs, support.f ~ -flange.f)
        compose(System(eqs, t; name = name), flange, support)
    else
        push!(eqs, s_support ~ 0)
        compose(System(eqs, t; name = name), flange)
    end
end

"""
    PartialElementaryTwoFlangesAndSupport2(; name, use_support = false)

Partial model for a component with two translational 1-dim. flanges and a support used for textual modeling, i.e., for elementary models

# Parameters:

  - `use_support`: If support flange enabled, otherwise implicitly grounded

# States:

  - `s_support`: [m] Absolute position of support flange"
"""
function PartialElementaryTwoFlangesAndSupport2(; name, use_support = false)
    @named flange = Flange()

    @variables s_a(t) [description = "Distance between left flange and support"]
    @variables s_b(t) [description = "Distance between right flange and support"]
    @variables s_support(t) [description = "Absolute position of support flange"]

    eqs = [s_a ~ flange_a.s - s_support
           s_b ~ flange_b.s - s_support]
    if use_support
        @named support = Support()
        push!(eqs, support.f ~ -flange_a.f - flange_b.f)
        compose(System(eqs, t; name = name), flange, support)
    else
        push!(eqs, s_support ~ 0)
        compose(System(eqs, t; name = name), flange)
    end
end

@mtkmodel PartialRigid begin
    @extend flange_a, flange_b = ptf = PartialTwoFlanges()
    @variables begin
        s(t), [description = "Absolute position of center of component"]
    end
    @parameters begin
        L = 0.0, [description = "Length of component, from left flange to right flange"]
    end
    @equations begin
        flange_a.s ~ s - L / 2
        flange_b.s ~ s + L / 2
    end
end
