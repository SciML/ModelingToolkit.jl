"""
    Ground(; name)

Ground node with the potential of zero and connector `g`. Every circuit must have one ground
node.

# Connectors:

  - `g`
"""
@mtkmodel Ground begin
    @components begin
        g = Pin()
    end
    @equations begin
        g.v ~ 0
    end
end

"""
    Resistor(; name, R = 1.0, T_ref = 300.15, alpha = 0, T_dep = false)

Generic resistor with optional temperature dependency.

# States:

  - See [OnePort](@ref)
  - `R(t)`: [`Ω`] Resistance (temperature dependent if `T_dep = true`)

# Connectors:

  - `p` Positive pin
  - `n` Negative pin
  - `heat_port` [HeatPort](@ref) (only if `T_dep = true`) Heat port to model the temperature dependency

# Parameters:

  - `R`: [`Ω`] Reference resistance
  - `T_ref`: [K] Reference temperature
  - `alpha`: [K⁻¹] Temperature coefficient of resistance
  - `T_dep`: [bool] Temperature dependency
"""
@mtkmodel Resistor begin
    @extend v, i = oneport = OnePort()

    @structural_parameters begin
        T_dep = false
    end

    @parameters begin
        R = 1.0, [description = "Reference resistance"]
        T_ref = 300.15, [description = "Reference temperature"]
        alpha = 0.0, [description = "Temperature coefficient of resistance"]
    end

    if T_dep
        @components begin
            heat_port = HeatPort()
        end
        @variables begin
            R_T(t), [description = "Temperature-dependent resistance"]
        end
        @equations begin
            R_T ~ R * (1 + alpha * (heat_port.T - T_ref))  # Temperature-dependent resistance
            heat_port.Q_flow ~ -v * i  # -LossPower
            v ~ i * R_T  # Ohm's Law
        end
    else
        @equations begin
            v ~ i * R  # Ohm's Law for constant resistance
        end
    end
end

"""
    Conductor(; name, G)

Creates an ideal conductor.

# States:

See [OnePort](@ref)

# Connectors:

  - `p` Positive pin
  - `n` Negative pin

# Parameters:

  - `G`: [`S`] Conductance
"""
@mtkmodel Conductor begin
    @extend v, i = oneport = OnePort()
    @parameters begin
        G, [description = "Conductance"]
    end
    @equations begin
        i ~ v * G
    end
end

"""
    Capacitor(; name, C, v)

Creates an ideal capacitor.
Initial voltage of capacitor can be set with `v` ([`V`])

# States:

See [OnePort](@ref)

# Connectors:

  - `p` Positive pin
  - `n` Negative pin

# Parameters:

  - `C`: [`F`] Capacitance
"""
@mtkmodel Capacitor begin
    @parameters begin
        C, [description = "Capacitance"]
    end
    @extend v, i = oneport = OnePort(; v)
    @equations begin
        D(v) ~ i / C
    end
end

"""
    Inductor(; name, L, i)

Creates an ideal Inductor.
Initial current through inductor can be set with `i` ([`A`]).

# States:

See [OnePort](@ref)

# Connectors:

  - `p` Positive pin
  - `n` Negative pin

# Parameters:

  - `L`: [`H`] Inductance
"""
@mtkmodel Inductor begin
    @parameters begin
        L, [description = "Inductance"]
    end
    @extend v, i = oneport = OnePort(; i)
    @equations begin
        D(i) ~ 1 / L * v
    end
end

"""
    IdealOpAmp(; name)

Ideal operational amplifier (norator-nullator pair).
The ideal OpAmp is a two-port. The left port is fixed to `v1 = 0` and `i1 = 0` (nullator).
At the right port both any voltage `v2` and any current `i2` are possible (norator).

# States:

See [TwoPort](@ref)

# Connectors:

  - `p1` Positive pin (left port)
  - `p2` Positive pin (right port)
  - `n1` Negative pin (left port)
  - `n2` Negative pin (right port)
"""
@mtkmodel IdealOpAmp begin
    @extend v1, v2, i1, i2 = twoport = TwoPort()
    @equations begin
        v1 ~ 0
        i1 ~ 0
    end
end

"""
    Short(; name)

Short is a simple short cut branch. That means the voltage drop between both pins is zero.

# States:

See [OnePort](@ref)

# Connectors:

  - `p` Positive pin
  - `n` Negative pin
"""
@mtkmodel Short begin
    @extend v, i = oneport = OnePort()
    @equations begin
        v ~ 0
    end
end

"""
    EMF(; name, k)

Electromotoric force (electric/mechanic transformer)

# States

  - `v(t)`: [`V`] The voltage across component `p.v - n.v`
  - `i(t)`: [`A`] The current passing through positive pin
  - `phi`: [`rad`] Rotation angle (=flange.phi - support.phi)
  - `w`: [`rad/s`] Angular velocity (= der(phi))

# Connectors

  - `p` [Pin](@ref) Positive pin
  - `n` [Pin](@ref) Negative pin
  - `flange` [Flange](@ref) Shaft of EMF shaft
  - `support` [Support](@ref) Support/housing of emf shaft

# Parameters:

  - `k`: [`N⋅m/A`] Transformation coefficient
"""
@mtkmodel EMF begin
    @parameters begin
        k, [description = "Transformation coefficient"]
    end
    @variables begin
        phi(t), [guess = 0.0, description = "Rotation Angle"]
        w(t), [guess = 0.0]
    end
    @extend v, i = oneport = OnePort()
    @components begin
        flange = Flange()
        support = Support()
    end
    @equations begin
        phi ~ flange.phi - support.phi
        D(phi) ~ w
        k * w ~ v
        flange.tau ~ -k * i
    end
end

"""
    Diode(; name, Is = 1e-6, n = 1, T = 300.15, T_dep = false)

Generic diode with optional temperature dependency.

# States

    - See [OnePort](@ref)

# Connectors

    - `p` Positive pin
    - `n` Negative pin
    - `port` [HeatPort](@ref) (only if `T_dep = true`) Heat port to model variable temperature dependency

# Parameters:
    
    - `Is`: [`A`] Saturation current
    - `n`: Ideality factor
    - `T`: [K] Constant ambient temperature - only used if T_dep=false
    - `T_dep`: [bool] Temperature dependency
"""
@mtkmodel Diode begin
    @constants begin
        k = 1.380649e-23 # Boltzmann constant (J/K)
        q = 1.602176634e-19 # Elementary charge (C)
    end
    @extend v, i = oneport = OnePort(; v = 0.0)

    @structural_parameters begin
        T_dep = false
    end

    @parameters begin
        Is = 1e-6, [description = "Saturation current (A)"]
        n = 1, [description = "Ideality factor"]
        T = 300.15, [description = "Ambient temperature"]
    end

    if T_dep
        @components begin
            port = HeatPort()
        end
        @variables begin
            Vt(t), [description = "Thermal voltage"]
        end
        @equations begin
            Vt ~ k * port.T / q  # Thermal voltage equation
            i ~ Is * (exp(v / (n * Vt)) - 1)  # Shockley diode equation with temperature dependence
            port.Q_flow ~ -v * i  # -LossPower
        end
    else
        @equations begin
            i ~ Is * (exp(v * q / (n * k * T)) - 1)  # Shockley diode equation
        end
    end
end

"""
    VariableResistor(; name, R_ref = 1.0, T_ref = 300.15, R_const = 1e-3, T_dep = false)

Variable resistor with optional temperature dependency.

The total resistance R ∈ [R_const, R_const + R_ref], where pos is the
position of the wiper and R_ref is the variable resistance between p and n.
The total resistance is then:

R = R_const + pos * R_ref

If T_dep is true, then R also depends on the temperature of the heat port with
temperature coefficient alpha. The total resistance is then:

R = R_const + pos * R_ref * (1 + alpha * (port.T - T_ref))

# States

    - See [OnePort](@ref)
    - `pos(t)`: Position of the wiper (normally 0-1)
    - `R(t)`: Resistance

# Connectors

        - `p` Positive pin
        - `n` Negative pin
        - `position` RealInput to set the position of the wiper
        - `port` [HeatPort](@ref) Heat port to model the temperature dependency

# Parameters

        - `R_ref`: [`Ω`] Resistance at temperature T_ref when fully closed (pos=1.0)
        - `T_ref`: [K] Reference temperature
        - `R_const`: [`Ω`] Constant resistance between p and n
        - `T_dep`: Temperature dependency
        - `alpha`: [K⁻¹] Temperature coefficient of resistance
        - `enforce_bounds`: Enforce bounds for the position of the wiper (0-1)
"""
@mtkmodel VariableResistor begin
    @extend v, i = oneport = OnePort()

    @structural_parameters begin
        T_dep = false
        enforce_bounds = true
    end

    @parameters begin
        R_ref = 1.0,
        [description = "Resistance at temperature T_ref when fully closed (pos=1.0) (Ω)"]
        T_ref = 300.15, [description = "Reference temperature (K)"]
        R_const = 1e-3, [description = "Constant resistance between p and n (Ω)"]
    end

    @components begin
        position = RealInput()
    end

    @variables begin
        pos(t), [description = "Position of the wiper (normally 0-1)"]
        R(t), [description = "Resistance (Ω)"]
    end

    if T_dep
        @parameters begin
            alpha = 1e-3, [description = "Temperature coefficient of resistance (K^-1)"]
        end
        @components begin
            port = HeatPort()
        end
        @equations begin
            port.Q_flow ~ -v * i  # -LossPower
            R ~ R_const + pos * R_ref * (1 + alpha * (port.T - T_ref))
        end
    else
        @equations begin
            R ~ R_const + pos * R_ref
        end
    end

    @equations begin
        pos ~ (enforce_bounds ? clamp(position.u, 0, 1) : position.u)
        v ~ i * R
    end
end
