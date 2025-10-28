"""
    MassFlow(; name)

Hydraulic mass flow input source

# Connectors:

  - `port`: hydraulic port
  - `dm`: real input
"""
@mtkmodel MassFlow begin
    @components begin
        port = HydraulicPort()
        dm = RealInput()
    end

    @equations begin
        port.dm ~ -dm.u
    end
end

"""
    FixedPressure(; p, name)

Fixed pressure source

# Parameters:
- `p`: [Pa] set pressure (set by `p` argument)

# Connectors:
- `port`: hydraulic port
"""
@mtkmodel FixedPressure begin
    @parameters begin
        p
    end

    @components begin
        port = HydraulicPort()
    end

    @equations begin
        port.p ~ p
    end
end
@deprecate Source FixedPressure

"""
    Pressure(; name)

input pressure source

# Connectors:
- `port`: hydraulic port
- `p`: real input
"""
@mtkmodel Pressure begin
    @components begin
        port = HydraulicPort()
        p = RealInput()
    end

    @equations begin
        port.p ~ p.u
    end
end
@deprecate InputSource Pressure
