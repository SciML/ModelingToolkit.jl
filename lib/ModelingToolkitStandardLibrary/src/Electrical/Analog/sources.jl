"""
    Voltage(; name)

Acts as an ideal voltage source with no internal resistance.

# States:

See [OnePort](@ref)

# Connectors:

  - `p` Positive pin
  - `n` Negative pin
  - `V` [RealInput](@ref) Input for the voltage control signal, i.e. `V ~ p.v - n.v`
"""
@mtkmodel Voltage begin
    @extend v, i = oneport = OnePort()
    @components begin
        V = RealInput()
    end
    @equations begin
        v ~ V.u
    end
end

"""
    Current(; name)

Acts as an ideal current source with no internal resistance.

# States:

See [OnePort](@ref)

# Connectors:

  - `p` Positive pin
  - `n` Negative pin
  - `I` [RealInput](@ref) Input for the current control signal, i.e. `I ~ p.i
"""
@mtkmodel Current begin
    @extend v, i = oneport = OnePort()
    @components begin
        I = RealInput()
    end
    @equations begin
        i ~ I.u
    end
end
