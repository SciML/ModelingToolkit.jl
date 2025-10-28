@connector MagneticPort begin
    V_m(t), [description = "Magnetic potential at the port"]
    Phi(t), [connect = Flow, description = "Magnetic flux flowing into the port"]
end
Base.@doc "Port for a Magnetic system." MagneticPort

"""
Positive magnetic port
"""
const PositiveMagneticPort = MagneticPort

"""
Negative magnetic port
"""
const NegativeMagneticPort = MagneticPort

"""
    TwoPort(; name, V_m = 0.0, Phi = 0.0)

Partial component with magnetic potential difference between two magnetic ports p and n and magnetic flux Phi from p to n.

# Parameters:

  - `V_m`: Initial magnetic potential difference between both ports
  - `Phi`: Initial magnetic flux from port_p to port_n
"""
@mtkmodel TwoPort begin
    @components begin
        port_p = PositiveMagneticPort()
        port_n = NegativeMagneticPort()
    end
    @variables begin
        V_m(t)
        Phi(t)
    end
    @equations begin
        V_m ~ port_p.V_m - port_n.V_m
        Phi ~ port_p.Phi
        0 ~ port_p.Phi + port_n.Phi
    end
end
