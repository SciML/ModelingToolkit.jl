"""
    AngleSensor(; name)

Ideal sensor to measure the absolute flange angle

# Connectors:

  - `flange`: [Flange](@ref) Flange of shaft from which sensor information shall be measured
  - `phi`: [RealOutput](@ref) Absolute angle of flange
"""
@mtkmodel AngleSensor begin
    @components begin
        flange = Flange()
        phi = RealOutput()
    end
    @equations begin
        phi.u ~ flange.phi
        flange.tau ~ 0
    end
end

"""
    SpeedSensor(; name)

Ideal sensor to measure the absolute flange angular velocity

# Connectors:

  - `flange`: [Flange](@ref) Flange of shaft from which sensor information shall be measured
  - `w`: [RealOutput](@ref) Absolute angular velocity of flange
"""
@mtkmodel SpeedSensor begin
    @components begin
        flange = Flange()
        w = RealOutput()
    end
    @equations begin
        D(flange.phi) ~ w.u
        flange.tau ~ 0
    end
end

"""
    TorqueSensor(;name)

Ideal sensor to measure the torque between two flanges (`= flange_a.tau`)

# Connectors:

  - `flange_a`: [Flange](@ref) Left flange of shaft
  - `flange_b`: [Flange](@ref) Left flange of shaft
  - `tau`: [RealOutput](@ref) Torque in flange flange_a and flange_b (`tau = flange_a.tau = -flange_b.tau`)
"""
@mtkmodel TorqueSensor begin
    @components begin
        flange_a = Flange()
        flange_b = Flange()
        tau = RealOutput()
    end
    @equations begin
        flange_a.phi ~ flange_b.phi
        tau.u ~ flange_a.tau
    end
end

"""
    RelSpeedSensor(; name)

Ideal sensor to measure the relative angular velocity

# Connectors:

  - `flange_a`: [Flange](@ref) Flange of shaft from which sensor information shall be measured
  - `flange_b`: [Flange](@ref) Flange of shaft from which sensor information shall be measured
  - `w`: [RealOutput](@ref) Absolute angular velocity of flange
"""
@mtkmodel RelSpeedSensor begin
    @components begin
        flange_a = Flange()
        flange_b = Flange()
        w_rel = RealOutput()
    end
    @variables begin
        phi_rel(t), [guess = 0.0]
    end
    @equations begin
        0 ~ flange_a.tau + flange_b.tau
        phi_rel ~ flange_b.phi - flange_a.phi
        D(phi_rel) ~ w_rel.u
        0 ~ flange_a.tau
    end
end
