"""
    ForceSensor(; name)

Linear 1D force sensor, measures the force between two flanges.

# Connectors:

- `flange_a`: 1-dim. translational flange
- `flange_b`: 1-dim. translational flange
- `output`: real output
"""
@mtkmodel ForceSensor begin
    @components begin
        flange_a = Flange()
        flange_b = Flange()
        output = RealOutput()
    end

    @equations begin
        flange_a.s ~ flange_b.s
        flange_a.f + flange_b.f ~ 0.0
        output.u ~ flange_a.f
    end
end

"""
    PositionSensor(; s = 0, name)

Linear 1D position sensor.

# States:

- `s`: [m] absolute position (with initial value of 0.0)

# Connectors:

- `flange`: 1-dim. translational flange
- `output`: real output
"""
@mtkmodel PositionSensor begin
    @components begin
        flange = Flange()
        output = RealOutput()
    end

    @equations begin
        output.u ~ flange.s
        flange.f ~ 0.0
    end
end

"""
    AccelerationSensor(; name)

Linear 1D acceleration sensor.

# States:

- `a`: [m/s^2] measured acceleration

# Connectors:

- `flange`: 1-dim. translational flange
- `output`: real output
"""
@mtkmodel AccelerationSensor begin
    @components begin
        flange = Flange()
        output = RealOutput()
    end

    @variables begin
        a(t)
    end

    @equations begin
        a ~ D(D(flange.s))
        output.u ~ a
        flange.f ~ 0.0
    end
end
