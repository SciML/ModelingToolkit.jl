@mtkmodel PartialTorque begin
    @extend flange,
    phi_support = partial_element = PartialElementaryOneFlangeAndSupport2(;
        use_support = false)
    @variables begin
        phi(t),
        [
            description = "Angle of flange with respect to support (= flange.phi - support.phi)"
        ]
    end
    @equations begin
        phi ~ flange.phi - phi_support
    end
end

"""
    Torque(; name, use_support = false)

Input signal acting as external torque on a flange

# States:

  - `phi_support(t)`: [`rad`] Absolute angle of support flange

# Connectors:

  - `flange` [Flange](@ref)
  - `tau` [RealInput](@ref)  Accelerating torque acting at flange `-flange.tau`

# Parameters:

  - `use_support`
"""
@mtkmodel Torque begin
    @extend (flange,) = partial_element = PartialElementaryOneFlangeAndSupport2(;
        use_support = false)
    @components begin
        tau = RealInput()
    end
    @equations begin
        flange.tau ~ -tau.u
    end
end

"""
    ConstantTorque(; name, tau_constant, use_support = false)

Constant torque source

# State variables:

- `phi_support(t)`: [`rad`] Absolute angle of support flange, only available if `use_support = true`
- `tau`: Accelerating torque acting at flange (= -flange.tau)
- `w`: Angular velocity of flange with respect to support (= der(phi))

# Connectors:
- `flange` [Flange](@ref)

# Arguments:
- `tau_constant`: The constant torque applied by the source
- `use_support`: Whether or not an internal support flange is added. By default, it is `false`
"""
@mtkmodel ConstantTorque begin
    @parameters begin
        tau_constant,
        [
            description = "Constant torque (if negative, torque is acting as load in positive direction of rotation)"
        ]
    end
    @extend flange, phi = partial_element = PartialTorque(; use_support = false)
    @variables begin
        tau(t), [description = "Accelerating torque acting at flange (= -flange.tau)"]
        w(t),
        [description = "Angular velocity of flange with respect to support (= der(phi))"]
    end
    @equations begin
        w ~ D(phi)
        tau ~ -flange.tau
        tau ~ tau_constant
    end
end

"""
    Speed(; name, use_support = false, exact = false, f_crit = 50)

Forced movement of a flange according to a reference angular velocity signal

# States:

  - `phi_support(t)`: [`rad`] Absolute angle of support flange"

# Connectors:

  - `flange` [Flange](@ref)
  - `w_ref` [RealInput](@ref) Reference angular velocity of flange with respect to support as input signal needs to be continuously differential

# Parameters:

  - `use_support`: If support flange enabled, otherwise implicitly grounded
  - `exact`: true/false exact treatment/filtering the input signal
  - `tau_filt`: [`rad/s`] if exact=false, Time constant of low-pass filter to filter input signal
"""
@component function Speed(; name, use_support = false, exact = false, tau_filt = 50)
    @named partial_element = PartialElementaryOneFlangeAndSupport2(use_support = use_support)
    @unpack flange, phi_support = partial_element
    @named w_ref = RealInput()
    @variables phi(t) [guess = 0.0] w(t) [guess = 0.0] a(t) [guess = 0.0]
    eqs = [phi ~ flange.phi - phi_support
           D(phi) ~ w]
    if exact
        pars = []
        push!(eqs, w ~ w_ref.u)
        push!(eqs, a ~ 0)
    else
        pars = @parameters tau_filt = tau_filt
        push!(eqs, D(w) ~ a)
        push!(eqs, a ~ (w_ref.u - w) * tau_filt)
    end
    return extend(System(eqs, t, [phi, w, a], pars; name = name, systems = [w_ref]),
        partial_element)
end

"""
    Position(; name, exact = false, f_crit = 50, use_support = false)

Forced movement of a flange according to a reference angle signal.

The input signal `phi_ref` defines the reference angle in [rad]. Flange is forced to move according to this reference motion relative to flange support. According to parameter `exact` (default = `false`), this is done in the following way:

- `exact=true`: The reference angle is treated exactly. This is only possible if the input signal is defined by an analytical function that can be differentiated at least twice in order to compute the acceleration.
- `exact=false`: The reference angle is filtered and the second derivative of the filtered curve is used to compute the reference acceleration of the flange. This second derivative is not computed by numerical differentiation but by an appropriate realization of the filter. For filtering, a second-order Bessel filter is used. The critical frequency (also called cut-off frequency) of the filter is defined via parameter `f_crit` in [Hz]. This value should be selected in such a way that it is higher than the essential low frequencies in the signal.

# Connectors 
- `flange::Flange`: Flange to be moved
- `phi_ref::RealInput`: Reference angle of flange with respect to support

# Variables
- `phi(t)`: Rotation angle of flange with respect to support
- `w(t)`: If `exact=false`, Angular velocity of flange with respect to support
- `a(t)`: If `exact=false`, Angular acceleration of flange with respect to support

# Parameters
- `exact`: (structural) true/false exact treatment/filtering the input signal
- `f_crit`: [Hz] if `exact=false`, Critical frequency of filter to filter input signal
"""
@component function Position(; name, exact = false, f_crit = 50, use_support = false)
    systems = @named begin
        partial_element = PartialElementaryOneFlangeAndSupport2(; use_support)
        phi_ref = RealInput()
    end
    @unpack flange, phi_support = partial_element

    pars = @parameters begin
        f_crit = f_crit, [description = "Critical frequency of input-signal filter"]
    end

    w_crit = 2 * π * f_crit
    af = 1.3617 # s coefficient of Bessel filter
    bf = 0.6180 # s*s coefficient of Bessel filter

    vars = @variables begin
        phi(t),
        [guess = 0.0, description = "Rotation angle of flange with respect to support"]
        w(t),
        [guess = 0.0, description = "Angular velocity of flange with respect to support"]
        a(t),
        [guess = 0.0,
            description = "Angular acceleration of flange with respect to support"]
    end

    equations = if exact
        [phi ~ flange.phi - phi_support
         phi ~ phi_ref.u]
    else
        [phi ~ flange.phi - phi_support
         D(phi) ~ w
         D(w) ~ a
         a ~ ((phi_ref.u - phi) * w_crit - af * w) * (w_crit / bf)]
    end
    extend(System(equations, t; name, systems = [phi_ref]), partial_element)
end
