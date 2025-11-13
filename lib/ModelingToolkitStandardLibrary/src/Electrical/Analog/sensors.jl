"""
    CurrentSensor(; name)

Creates a circuit component that measures the current flowing through it. Analogous to
an ideal ammeter.

# States:

  - `i(t)`: [`A`] Current through the sensor

# Connectors:

  - `p` Positive pin
  - `n` Negative pin
"""
@mtkmodel CurrentSensor begin
    @components begin
        p = Pin()
        n = Pin()
    end
    @variables begin
        i(t)
    end
    @equations begin
        p.v ~ n.v
        i ~ p.i
        i ~ -n.i
    end
end

"""
PotentialSensor(; name)

Creates a circuit component which measures the potential at a pin.

# States:

  - `phi(t)`: [`V`] The measured potential at this point

# Connectors:

  - `p` Pin at which potential is to be measured
"""
@mtkmodel PotentialSensor begin
    @components begin
        p = Pin()
    end
    @variables begin
        phi(t)
    end
    @equations begin
        p.i ~ 0
        phi ~ p.v
    end
end

"""
VoltageSensor(; name)

Creates a circuit component that measures the voltage across it. Analogous to an ideal voltmeter.

# States:

  - `v(t)`: [`V`] The voltage difference from positive to negative pin `p.v - n.v`

# Connectors:

  - `p` Positive pin
  - `n` Negative pin
"""
@mtkmodel VoltageSensor begin
    @components begin
        p = Pin()
        n = Pin()
    end
    @variables begin
        v(t)
    end
    @equations begin
        p.i ~ 0
        n.i ~ 0
        v ~ p.v - n.v
    end
end

"""
PowerSensor(; name)

Combines a [`VoltageSensor`](@ref) and a [`CurrentSensor`](@ref) to measure the power being
consumed by a circuit.

# States:

  - `power(t)`: [`W`] The power being consumed, given by the product of voltage and current
  - See [VoltageSensor](@ref)
  - See [CurrentSensor](@ref)

# Connectors:

  - `pc` Corresponds to the `p` pin of the [`CurrentSensor`](@ref)
  - `nc` Corresponds to the `n` pin of the [`CurrentSensor`](@ref)
  - `pv` Corresponds to the `p` pin of the [`VoltageSensor`](@ref)
  - `nv` Corresponds to the `n` pin of the [`VoltageSensor`](@ref)
"""
@mtkmodel PowerSensor begin
    @components begin
        pc = Pin()
        nc = Pin()
        pv = Pin()
        nv = Pin()
        voltage_sensor = VoltageSensor()
        current_sensor = CurrentSensor()
    end
    @variables begin
        power(t)
    end
    @equations begin
        connect(voltage_sensor.p, pv)
        connect(voltage_sensor.n, nv)
        connect(current_sensor.p, pc)
        connect(current_sensor.n, nc)
        power ~ current_sensor.i * voltage_sensor.v
    end
end

"""
MultiSensor(; name)

Combines a [`VoltageSensor`](@ref) and a [`CurrentSensor`](@ref).

# States:

  - `v(t)`: [`V`] The voltage across the [`VoltageSensor`](@ref). Defaults to 1.0.
  - `i(t)`: [`A`] The current across the [`CurrentSensor`](@ref). Defaults to 1.0.

# Connectors:

  - `pc` Corresponds to the `p` pin of the [`CurrentSensor`](@ref)
  - `nc` Corresponds to the `n` pin of the [`CurrentSensor`](@ref)
  - `pv` Corresponds to the `p` pin of the [`VoltageSensor`](@ref)
  - `nv` Corresponds to the `n` pin of the [`VoltageSensor`](@ref)
"""
@mtkmodel MultiSensor begin
    @components begin
        pc = Pin()
        nc = Pin()
        pv = Pin()
        nv = Pin()
        voltage_sensor = VoltageSensor()
        current_sensor = CurrentSensor()
    end
    @variables begin
        i(t) = 1.0
        v(t) = 1.0
    end
    @equations begin
        connect(voltage_sensor.p, pv)
        connect(voltage_sensor.n, nv)
        connect(current_sensor.p, pc)
        connect(current_sensor.n, nc)
        i ~ current_sensor.i
        v ~ voltage_sensor.v
    end
end
