"""
    Fixed(; name, s0 = 0.0)

Flange fixed in housing at a given position.

# Parameters:

  - `s0`: [m] Fixed offset position of housing

# Connectors:

  - `flange: 1-dim. translational flange`
"""
@mtkmodel Fixed begin
    @parameters begin
        s0 = 0
    end

    @components begin
        flange = Flange()
    end

    @equations begin
        flange.s ~ s0
    end
end

"""
    Mass(; name, m, s, v = 0.0)

Sliding mass with inertia

# Parameters:

  - `m`: [kg] Mass of sliding mass

# States:

  - `s`: [m] Absolute position of sliding mass. It accepts an initial value, which defaults to 0.0.
  - `v`: [m/s] Absolute linear velocity of sliding mass (= D(s)). It accepts an initial value, which defaults to 0.0.

# Connectors:

  - `flange: 1-dim. translational flange of mass`
"""
@mtkmodel Mass begin
    @parameters begin
        m = 0.0, [description = "Mass of sliding mass [kg]"]
    end
    @variables begin
        v(t), [description = "Absolute linear velocity of sliding mass [m/s]"]
        a(t), [description = "Absolute linear acceleration of sliding mass [m/s^2]"]
    end
    @extend flange_a, flange_b, s = pr = PartialRigid(; L = 0.0, s)
    @equations begin
        v ~ D(s)
        a ~ D(v)
        m * a ~ flange_a.f + flange_b.f
    end
end

"""
    Spring(; c= 0.0, name, s_rel0 = 0)

Linear 1D translational spring

# Parameters:

  - `c`: [N/m] Spring constant
  - `s_rel0`: Unstretched spring length

# Connectors:

  - `flange_a: 1-dim. translational flange on one side of spring`
  - `flange_b: 1-dim. translational flange on opposite side of spring` #default function
"""
@mtkmodel Spring begin
    @extend flange_a, flange_b, s_rel, f = pc = PartialCompliant()
    @parameters begin
        c = 0.0, [description = "Spring constant [N/m]"]
        s_rel0 = 0.0, [description = "Unstretched spring length [m]"]
    end

    @equations begin
        f ~ c * (s_rel - s_rel0)
    end
end

"""
    Damper(; name, d = 0.0)

Linear 1D translational damper

# Parameters:

  - `d`: [N.s/m] Damping constant

# Connectors:

  - `flange_a: 1-dim. translational flange on one side of damper`
  - `flange_b: 1-dim. translational flange on opposite side of damper`
"""
@mtkmodel Damper begin
    @extend flange_a, flange_b, v_rel, f = pc = PartialCompliantWithRelativeStates()
    @parameters begin
        d = 0.0, [description = "Damping constant [Ns/m]"]
    end
    @variables begin
        lossPower(t), [description = "Power dissipated by the damper [W]"]
    end
    @equations begin
        f ~ d * v_rel
        lossPower ~ f * v_rel
    end
end

"""
    SpringDamper(; name, c = 0.0, d = 0.0, s_rel0 = 0.0)

Linear 1D translational spring and damper in parallel

# Parameters:
- `c`: [N/m] Spring constant
- `d`: [N.s/m] Damping constant
- `s_rel0`: Unstretched spring length

# Connectors:
- `flange_a: 1-dim. translational flange on one side of spring`
- `flange_b: 1-dim. translational flange on opposite side of spring`

# Variables:
- `lossPower`: [W] Power dissipated by the damper
- `f`: [N] Total force
"""
@mtkmodel SpringDamper begin
    @extend flange_a, flange_b, s_rel, v_rel, f = pc = PartialCompliantWithRelativeStates()
    @parameters begin
        d = 0.0, [description = "Damping constant [Ns/m]"]
        c = 0.0, [description = "Spring constant [N/m]"]
        s_rel0 = 0.0, [description = "Unstretched spring length [m]"]
    end
    @variables begin
        lossPower(t), [description = "Power dissipated by the damper [W]"]
    end
    @equations begin
        f ~ c * (s_rel - s_rel0) + d * v_rel
        lossPower ~ d * v_rel^2
    end
end
