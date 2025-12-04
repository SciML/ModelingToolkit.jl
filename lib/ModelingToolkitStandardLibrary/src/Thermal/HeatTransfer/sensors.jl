"""
    TemperatureSensor(; name)

Absolute temperature sensor in kelvin.

This is an ideal absolute temperature sensor which returns the temperature of the connected port in kelvin as an output
signal. The sensor itself has no thermal interaction with whatever it is connected to. Furthermore, no thermocouple-like
lags are associated with this sensor model.

# Connectors:

  - `port`: [HeatPort](@ref) Thermal port from which sensor information shall be measured
  - `T`: [RealOutput](@ref) [K] Absolute temperature of port
"""
@mtkmodel TemperatureSensor begin
    @components begin
        port = HeatPort()
        T = RealOutput()
    end
    @equations begin
        T.u ~ port.T
        port.Q_flow ~ 0
    end
end

"""
    RelativeTemperatureSensor(; name)

Relative Temperature sensor.

The relative temperature `port_a.T - port_b.T` is determined between the two ports of this component and is provided as
output signal in kelvin.

# Connectors:

  - `port_a`: [HeatPort](@ref) Thermal port from which sensor information shall be measured
  - `port_b`: [HeatPort](@ref) Thermal port from which sensor information shall be measured
  - `T`: [RealOutput](@ref) [K] Relative temperature `a.T - b.T`
"""
@mtkmodel RelativeTemperatureSensor begin
    @components begin
        port_a = HeatPort()
        port_b = HeatPort()
        T = RealOutput()
    end
    @equations begin
        T.u ~ port_a.T - port_b.T
        port_a.Q_flow ~ 0
        port_b.Q_flow ~ 0
    end
end

"""
    HeatFlowSensor(; name)

Heat flow rate sensor.

This model is capable of monitoring the heat flow rate flowing through this component. The sensed value of heat flow rate
is the amount that passes through this sensor while keeping the temperature drop across the sensor zero. This is an ideal
model, so it does not absorb any energy, and it has no direct effect on the thermal response of a system it is included in.
The output signal is positive, if the heat flows from `port_a` to `port_b`.

# Connectors:

  - `port_a`: [HeatPort](@ref) Thermal port from which sensor information shall be measured
  - `port_b`: [HeatPort](@ref) Thermal port from which sensor information shall be measured
  - `Q_flow`: [RealOutput](@ref) [W] Heat flow from `port_a` to `port_b` 
"""
@mtkmodel HeatFlowSensor begin
    @components begin
        port_a = HeatPort()
        port_b = HeatPort()
        Q_flow = RealOutput()
    end
    @equations begin
        port_a.T ~ port_b.T
        port_a.Q_flow + port_b.Q_flow ~ 0
        Q_flow.u ~ port_a.Q_flow
    end
end
