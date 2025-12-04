"""
    NMOS(;name, V_tn, R_DS, lambda)

Creates an N-type MOSFET transistor 

    # Structural Parameters
        - `use_transconductance`: If `true` the parameter `k_n` needs to be provided, and is used in the calculation of the current
        through the transistor. Otherwise, `mu_n`, `C_ox`, `W`, and `L` need to be provided and are used to calculate the transconductance. 

        - `use_channel_length_modulation`: If `true` the channel length modulation effect is taken in to account. In essence this gives  
        the drain-source current has a small dependency on the drains-source voltage in the saturation region of operation. 

    # Connectors
        - `d` Drain Pin
        - `g` Gate Pin
        - `s` Source Pin

    # Parameters
        - `mu_n`: Electron mobility
        - `C_ox`: Oxide capacitance (F/m^2)
        - `W`: Channel width (m)
        - `L`: Channel length
        - `k_n`: MOSFET transconductance parameter 

Based on the MOSFET models in (Sedra, A. S., Smith, K. C., Carusone, T. C., & Gaudet, V. C. (2021). Microelectronic circuits (8th ed.). Oxford University Press.)
"""
@mtkmodel NMOS begin
    @variables begin
        V_GS(t)
        V_DS(t)
        V_OV(t)
    end

    @components begin
        d = Pin()
        g = Pin()
        s = Pin()
    end

    @parameters begin
        V_tn = 0.8, [description = "Threshold voltage (V)"]
        R_DS = 1e7, [description = "Drain to source resistance (Ω)"]

        lambda = 0.04, [description = "Channel length modulation coefficient (V^(-1))"]

        if !use_transconductance
            mu_n, [description = "Electron mobility"]
            C_ox, [description = "Oxide capacitance (F/m^2)"]
            W, [description = "Channel width (m)"]
            L, [description = "Channel length (m)"]
        else
            k_n = 20e-3, [description = "MOSFET transconductance parameter"]
        end
    end

    @structural_parameters begin
        use_transconductance = true
    end

    begin
        if !use_transconductance
            k_n = mu_n * C_ox * (W / L)
        end
    end

    @equations begin
        V_DS ~ ifelse(d.v < s.v, s.v - d.v, d.v - s.v)
        V_GS ~ g.v - ifelse(d.v < s.v, d.v, s.v)
        V_OV ~ V_GS - V_tn

        d.i ~
        ifelse(d.v < s.v, -1, 1) * ifelse(V_GS < V_tn,
            V_DS / R_DS,
            ifelse(V_DS < V_OV,
                k_n * (1 + lambda * V_DS) * (V_OV - V_DS / 2) * V_DS + V_DS / R_DS,
                ((k_n * V_OV^2) / 2) * (1 + lambda * V_DS) + V_DS / R_DS
            )
        )

        g.i ~ 0
        s.i ~ -d.i
    end
end

"""
    PMOS(;name, V_tp, R_DS, lambda)

Creates an N-type MOSFET transistor 

    # Structural Parameters
        - `use_transconductance`: If `true` the parameter `k_p` needs to be provided, and is used in the calculation of the current
        through the transistor. Otherwise, `mu_n`, `C_ox`, `W`, and `L` need to be provided and are used to calculate the transconductance. 

        - `use_channel_length_modulation`: If `true` the channel length modulation effect is taken in to account. In essence this gives  
        the drain-source current has a small dependency on the drains-source voltage in the saturation region of operation. 

    # Connectors
        - `d` Drain Pin
        - `g` Gate Pin
        - `s` Source Pin

    # Parameters
        - `mu_p`: Electron mobility
        - `C_ox`: Oxide capacitance (F/m^2)
        - `W`: Channel width (m)
        - `L`: Channel length
        - `k_p`: MOSFET transconductance parameter 

Based on the MOSFET models in (Sedra, A. S., Smith, K. C., Carusone, T. C., & Gaudet, V. C. (2021). Microelectronic circuits (8th ed.). Oxford University Press.)
"""
@mtkmodel PMOS begin
    @variables begin
        V_GS(t)
        V_DS(t)
    end

    @components begin
        d = Pin()
        g = Pin()
        s = Pin()
    end

    @parameters begin
        V_tp = -1.5, [description = "Threshold voltage (V)"]
        R_DS = 1e7, [description = "Drain-source resistance (Ω)"]

        lambda = 1 / 25, [description = "Channel length modulation coefficient (V^(-1))"]

        if !use_transconductance
            mu_p, [description = "Hole mobility"]
            C_ox, [description = "Oxide capacitance (F/m^2)"]
            W, [description = "Channel width (m)"]
            L, [description = "Channel length (m)"]
        else
            k_p = 20e-3
        end
    end

    @structural_parameters begin
        use_transconductance = true
    end

    begin
        if !use_transconductance
            k_p = mu_p * C_ox * (W / L)
        end
    end

    @equations begin
        V_DS ~ ifelse(d.v > s.v, s.v - d.v, d.v - s.v)
        V_GS ~ g.v - ifelse(d.v > s.v, d.v, s.v)

        d.i ~
        -ifelse(d.v > s.v, -1.0, 1.0) * ifelse(V_GS > V_tp,
            V_DS / R_DS,
            ifelse(V_DS > (V_GS - V_tp),
                k_p * (1 + lambda * V_DS) * ((V_GS - V_tp) - V_DS / 2) * V_DS +
                V_DS / R_DS,
                ((k_p * (V_GS - V_tp)^2) / 2) * (1 + lambda * V_DS) + V_DS / R_DS
            )
        )

        g.i ~ 0
        s.i ~ -d.i
    end
end
