"""
    ConstantMagneticPotentialDifference(; name, V_m = 0.0)

Constant magnetomotive force.

Parameters:

  - `V_m`: [A] Magnetic potential difference
"""
@mtkmodel ConstantMagneticPotentialDifference begin
    @components begin
        port_p = PositiveMagneticPort()
        port_n = NegativeMagneticPort()
    end
    @parameters begin
        V_m = 0.0, [description = "Magnetic potential difference"]
    end
    @variables begin
        Phi(t)
    end
    @equations begin
        V_m ~ port_p.V_m - port_n.V_m
        Phi ~ port_p.Phi
        0 ~ port_p.Phi + port_n.Phi
    end
end

"""
    ConstantMagneticFlux(; name, Phi = 0.0)

Source of constant magnetic flux.

Parameters:

  - `Phi`: [Wb] Magnetic flux
"""
@mtkmodel ConstantMagneticFlux begin
    @components begin
        port_p = PositiveMagneticPort()
        port_n = NegativeMagneticPort()
    end
    @parameters begin
        Phi = 0.0, [description = "Magnetic flux"]
    end
    @variables begin
        V_m(t)
    end
    @equations begin
        V_m ~ port_p.V_m - port_n.V_m
        Phi ~ port_p.Phi
        0 ~ port_p.Phi + port_n.Phi
    end
end
