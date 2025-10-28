"""
    Fixed(;name, phi0 = 0.0)

Flange fixed in housing at a given angle.

# Connectors:

  - `flange` [Flange](@ref)

# Parameters:

  - `phi0`: [`rad`] Fixed offset angle of housing
"""
@mtkmodel Fixed begin
    @components begin
        flange = Flange()
    end
    @parameters begin
        phi0 = 0.0, [description = "Fixed offset angle of flange"]
    end
    @equations begin
        flange.phi ~ phi0
    end
end

"""
    Inertia(;name, J, phi = 0.0, w = 0.0, a = 0.0)

1D-rotational component with inertia.

# States:

  - `phi`: [`rad`] Absolute rotation angle of component
  - `w`: [`rad/s`] Absolute angular velocity of component (= D(phi))
  - `a`: [`rad/s²`] Absolute angular acceleration of component (= D(w))

# Connectors:

  - `flange_a` [Flange](@ref) Left flange
  - `flange_b` [Flange](@ref) Right flange

# Parameters:

  - `J`: [`kg·m²`] Moment of inertia
"""
@mtkmodel Inertia begin
    @parameters begin
        J, [description = "Moment of inertia"]
    end
    @components begin
        flange_a = Flange()
        flange_b = Flange()
    end
    begin
        @symcheck J > 0 || throw(ArgumentError("Expected `J` to be positive"))
    end
    @variables begin
        phi(t), [description = "Absolute rotation angle", guess = 0.0]
        w(t), [description = "Absolute angular velocity", guess = 0.0]
        a(t), [description = "Absolute angular acceleration", guess = 0.0]
    end
    @equations begin
        phi ~ flange_a.phi
        phi ~ flange_b.phi
        D(phi) ~ w
        D(w) ~ a
        J * a ~ flange_a.tau + flange_b.tau
    end
end

"""
    Spring(; name, c, phi_rel0 = 0.0)

Linear 1D rotational spring

# States:

  - `phi_rel(t)`: [`rad`] Relative rotation angle (`flange_b.phi - flange_a.phi`)
  - `tau(t)`: [`N.m`] Torque between flanges (`flange_b.tau`)

# Connectors:

  - `flange_a` [Flange](@ref)
  - `flange_b` [Flange](@ref)

# Parameters:

  - `c`: [`N.m/rad`] Spring constant
  - `phi_rel0`: [`rad`] Unstretched spring angle. Defaults to 0.0.
"""
@mtkmodel Spring begin
    @extend phi_rel, tau = partial_comp = PartialCompliant()
    begin
        @symcheck c > 0 || throw(ArgumentError("Expected `c` to be positive"))
    end
    @parameters begin
        c, [description = "Spring constant"]
        phi_rel0 = 0.0, [description = "Unstretched spring angle"]
    end
    @equations begin
        tau ~ c * (phi_rel - phi_rel0)
    end
end

"""
    Damper(; name, d)

Linear 1D rotational damper

# States:

  - `phi_rel(t)`: [`rad`] Relative rotation angle (= flange_b.phi - flange_a.phi)
  - `w_rel(t)`: [`rad/s`] Relative angular velocity (= D(phi_rel))
  - `a_rel(t)`: [`rad/s²`] Relative angular acceleration (= D(w_rel))
  - `tau(t)`: [`N.m`] Torque between flanges (= flange_b.tau)

# Connectors:

  - `flange_a` [Flange](@ref)
  - `flange_b` [Flange](@ref)

# Parameters:

  - `d`: [`N.m.s/rad`] Damping constant
"""
@mtkmodel Damper begin
    @extend w_rel, tau = partial_comp = PartialCompliantWithRelativeStates()
    begin
        @symcheck d > 0 || throw(ArgumentError("Expected `d` to be positive"))
    end
    @parameters begin
        d, [description = "Damping constant"]
    end
    @equations begin
        tau ~ d * w_rel
    end
end
"""
    SpringDamper(; name, d)

Linear 1D rotational spring and damper

# States:

  - `phi_rel(t)`: [`rad`] Relative rotation angle (= flange_b.phi - flange_a.phi)
  - `w_rel(t)`: [`rad/s`] Relative angular velocity (= D(phi_rel))
  - `a_rel(t)`: [`rad/s²`] Relative angular acceleration (= D(w_rel))
  - `tau(t)`: [`N.m`] Torque between flanges (= flange_b.tau)

# Connectors:

  - `flange_a` [Flange](@ref)
  - `flange_b` [Flange](@ref)

# Parameters:

  - `d`: [`N.m.s/rad`] Damping constant
  - `c`: [`N.m/rad`] Spring constant
  - `phi_rel0`: [`rad`] Unstretched spring angle. Defaults to 0.0
"""
@mtkmodel SpringDamper begin
    @extend phi_rel, w_rel, tau = partial_comp = PartialCompliantWithRelativeStates()
    @variables begin
        tau_c(t), [description = "Spring torque"]
        tau_d(t), [description = "Damper torque"]
    end
    @parameters begin
        d, [description = "Damping constant"]
        c, [description = "Spring constant"]
        phi_rel0 = 0.0, [description = "Unstretched spring angle"]
    end
    @equations begin
        tau_c ~ c * (phi_rel - phi_rel0)
        tau_d ~ d * w_rel
        tau ~ tau_c + tau_d
    end
end

"""
    IdealGear(; name, ratio, use_support = false)

Ideal gear without inertia.

This element characterizes any type of gear box which is fixed in the ground and which has one driving shaft and one driven shaft.

# States:

  - `phi_a(t)`: [`rad`] Relative angle between shaft a and the support
  - `phi_b(t)`: [`rad`] Relative angle between shaft b and the support

# Connectors:

  - `flange_a` [Flange](@ref)
  - `flange_b` [Flange](@ref)
  - `support` [Support](@ref) if `use_support == true`

# Parameters:

  - `ratio`: Transmission ratio (flange_a.phi/flange_b.phi)
  - `use_support`: If support flange enabled, otherwise implicitly grounded. By default it is `false`
"""
@mtkmodel IdealGear begin
    @extend phi_support, flange_a,
    flange_b = partial_element = PartialElementaryTwoFlangesAndSupport2(;
        use_support = false)
    @parameters begin
        ratio, [description = "Transmission ratio"]
    end
    @variables begin
        phi_a(t),
        [description = "Relative angle between shaft a and the support", guess = 0.0]
        phi_b(t),
        [description = "Relative angle between shaft b and the support", guess = 0.0]
    end
    @equations begin
        phi_a ~ flange_a.phi - phi_support
        phi_b ~ flange_b.phi - phi_support
        phi_a ~ ratio * phi_b
        0 ~ ratio * flange_a.tau + flange_b.tau
    end
end

"""
    RotationalFriction(; name, f, tau_c, w_brk, tau_brk)

Models rotational friction with Stribeck effect, Coulomb friction and viscous friction between the two flanges.
The friction torque is a function of the relative angular velocity between `flange_a` and `flange_b`.

Friction model: "Armstrong, B. and C.C. de Wit, Friction Modeling and Compensation, The Control Handbook, CRC Press, 1995."

# States:

  - `phi_rel(t)`: [`rad`] Relative rotation angle `(= flange_b.phi - flange_a.phi)`
  - `w_rel(t)`: [`rad/s`] Relative angular velocity `(= D(phi_rel))`
  - `a_rel(t)`: [`rad/s²`] Relative angular acceleration `(= D(w_rel))`
  - `tau(t)`: [`N.m`] Torque between flanges `(= flange_b.tau)`

# Connectors:

  - `flange_a` [Flange](@ref)
  - `flange_b` [Flange](@ref)

# Parameters:

  - `f`: [`N⋅m/(rad/s)`] Viscous friction coefficient
  - `tau_c`: [`N⋅m`] Coulomb friction torque
  - `w_brk`: [`rad/s`] Breakaway friction velocity
  - `tau_brk`: [`N⋅m`] Breakaway friction torque
"""
@mtkmodel RotationalFriction begin
    @extend w_rel, tau = partial_comp = PartialCompliantWithRelativeStates()
    @parameters begin
        f, [description = "Viscous friction coefficient"]
        tau_c, [description = "Coulomb friction torque"]
        w_brk, [description = "Breakaway friction velocity"]
        tau_brk, [description = "Breakaway friction torque"]
    end

    begin
        str_scale = sqrt(2 * exp(1)) * (tau_brk - tau_c)
        w_st = w_brk * sqrt(2)
        w_coul = w_brk / 10
    end
    @equations begin
        tau ~
        str_scale * (exp(-(w_rel / w_st)^2) * w_rel / w_st) +
        tau_c * tanh(w_rel / w_coul) + f * w_rel # Stribeck friction + Coulomb friction + Viscous friction
    end
end
