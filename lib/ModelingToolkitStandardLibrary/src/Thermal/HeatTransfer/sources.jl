"""
    FixedHeatFlow(; name, Q_flow = 1.0, T_ref = 293.15, alpha = 0.0)

Fixed heat flow boundary condition.

This model allows a specified amount of heat flow rate to be "injected" into a thermal system at a given port.
The constant amount of heat flow rate `Q_flow` is given as a parameter. The heat flows into the component to which
the component FixedHeatFlow is connected, if parameter `Q_flow` is positive.

# Connectors:

  - `port`

# Parameters:

  - `Q_flow`: [W] Fixed heat flow rate at port
  - `T_ref`: [K] Reference temperature
  - `alpha`: [1/K] Temperature coefficient of heat flow rate
"""
@mtkmodel FixedHeatFlow begin
    @parameters begin
        Q_flow = 1.0, [description = "Fixed heat flow rate at port"]
        T_ref = 293.15, [description = "Reference temperature"]
        alpha = 0.0, [description = "Temperature coefficient of heat flow rate"]
    end
    @components begin
        port = HeatPort()
    end

    @equations begin
        port.Q_flow ~ ifelse(alpha == 0.0,
            -Q_flow,  # Simplified equation when alpha is 0
            -Q_flow * (1 + alpha * (port.T - T_ref)))
    end
end

"""
    FixedTemperature(; name, T)

Fixed temperature boundary condition in kelvin.

This model defines a fixed temperature `T` at its port in kelvin, i.e., it defines a fixed temperature as a boundary condition.

# Connectors:

  - `port`

# Parameters:

  - `T`: [K] Fixed temperature boundary condition
"""
@mtkmodel FixedTemperature begin
    @components begin
        port = HeatPort()
    end
    @parameters begin
        T, [description = "Fixed temperature boundary condition"]
    end
    @equations begin
        port.T ~ T
    end
end

"""
    PrescribedHeatFlow(; name, T_ref = 293.15, alpha = 0.0)

Prescribed heat flow boundary condition.

This model allows a specified amount of heat flow rate to be "injected" into a thermal system at a given port.
The amount of heat is given by the input signal `Q_flow` into the model. The heat flows into the component to which
the component `PrescribedHeatFlow` is connected, if the input signal is positive.
If parameter alpha is > 0, the heat flow is multiplied by `1 + alpha*(port.T - T_ref`) in order to simulate temperature
dependent losses (which are given a reference temperature T_ref).

# Connectors:

  - `port`
  - `RealInput` `Q_flow` Input for the heat flow

# Parameters:

  - `T_ref`: [K] Reference temperature
  - `alpha`: [1/K] Temperature coefficient of heat flow rate
"""
@mtkmodel PrescribedHeatFlow begin
    @parameters begin
        T_ref = 293.15, [description = "Reference temperature"]
        alpha = 0.0, [description = "Temperature coefficient of heat flow rate"]
    end
    @components begin
        port = HeatPort()
        Q_flow = RealInput()
    end
    @equations begin
        port.Q_flow ~ -Q_flow.u * (1 + alpha * (port.T - T_ref))
    end
end

"""
    PrescribedTemperature(; name)

This model represents a variable temperature boundary condition.

The temperature in kelvin is given as input signal to the `RealInput` `T`. The effect is that an instance of
this model acts as an infinite reservoir, able to absorb or generate as much energy as required to keep
the temperature at the specified value.

# Connectors:

  - `port`
  - `RealInput` `T` input for the temperature
"""
@mtkmodel PrescribedTemperature begin
    @components begin
        port = HeatPort()
        T = RealInput()
    end
    @equations begin
        port.T ~ T.u
    end
end
