"""
    Ground(; name)

Zero magnetic potential.
"""
@mtkmodel Ground begin
    @components begin
        port = PositiveMagneticPort()
    end
    @equations begin
        port.V_m ~ 0
    end
end

"""
    Idle(;name)

Idle running branch.
"""
@mtkmodel Idle begin
    @extend (Phi,) = two_port = TwoPort()
    @equations begin
        Phi ~ 0
    end
end

"""
    Short(;name)

Short cut branch.
"""
@mtkmodel Short begin
    @extend (V_m,) = two_port = TwoPort()
    @equations begin
        V_m ~ 0
    end
end

"""
    Crossing(;name)

Crossing of two branches.

This is a simple crossing of two branches. The ports port_p1 and port_p2 are connected, as well as port_n1 and port_n2.
"""
@mtkmodel Crossing begin
    @components begin
        port_p1 = PositiveMagneticPort()
        port_p2 = PositiveMagneticPort()
        port_n1 = NegativeMagneticPort()
        port_n2 = NegativeMagneticPort()
    end
    @equations begin
        connect(port_p1, port_p2)
        connect(port_n1, port_n2)
    end
end

"""
    ConstantPermeance(; name, G_m = 1.0)

Constant permeance.

# Parameters:

  - `G_m`: [H] Magnetic permeance
"""
@mtkmodel ConstantPermeance begin
    @extend V_m, Phi = two_port = TwoPort()
    @parameters begin
        G_m = 1.0, [description = "Magnetic permeance"]
    end
    @equations begin
        Phi ~ G_m * V_m
    end
end

"""
    ConstantReluctance(; name, R_m = 1.0)

Constant reluctance.

# Parameters:

  - `R_m`: [H^-1] Magnetic reluctance
"""
@mtkmodel ConstantReluctance begin
    @extend V_m, Phi = two_port = TwoPort(; Phi = 0.0)
    @parameters begin
        R_m = 1.0, [description = "Magnetic reluctance"]
    end
    @equations begin
        V_m ~ Phi * R_m
    end
end

"""
    ElectroMagneticConverter(; name, N, Phi)

Ideal electromagnetic energy conversion.

The electromagnetic energy conversion is given by Ampere's law and Faraday's law respectively
V_m = N * i
N * dΦ/dt = -v

Initial magnetic flux flowing into the port_p can be set with `Phi` ([Wb])

# Parameters:

  - `N`: Number of turns
"""
@mtkmodel ElectroMagneticConverter begin
    @parameters begin
        N, [description = "Number of turns"]
    end
    @variables begin
        v(t)
        i(t)
    end
    @extend V_m, Phi = two_port = TwoPort(; Phi)
    @components begin
        p = Pin()
        n = Pin()
    end
    @equations begin
        v ~ p.v - n.v
        0 ~ p.i + n.i
        i ~ p.i
        #converter equations:
        V_m ~ i * N # Ampere's law
        D(Phi) ~ -v / N
    end
end

"""
    EddyCurrent(;name, Phi, rho = 0.098e-6, l = 1, A = 1)

For modelling of eddy current in a conductive magnetic flux tube.
Initial magnetic flux flowing into the port_p can be set with `Phi` ([`Wb`])

# Parameters:

  - `rho`: [ohm * m] Resistivity of flux tube material (default: Iron at 20degC)
  - `l`: [m] Average length of eddy current path
  - `A`: [m^2] Cross sectional area of eddy current path
"""
@mtkmodel EddyCurrent begin
    @parameters begin
        rho = 0.098e-6, [description = "Resistivity of flux tube material"]
        l = 1, [description = "Average length of eddy current path"]
        A = 1, [description = "Cross sectional area of eddy current path"]
        R = rho * l / A # Electrical resistance of eddy current path
    end
    @extend (V_m, Phi) = two_port = TwoPort(; Phi)
    @equations begin
        D(Phi) ~ V_m * R
    end
end
