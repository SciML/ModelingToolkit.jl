"""
    NPN(;name, B_F, B_R, Is, V_T, V_A, phi_C, phi_E, Z_C, Z_E, Tau_f, Tau_r, C_jC0, C_jE0, C_CS, gamma_C, gamma_E, NF, NR)

Creates an NPN Bipolar Junction Transistor following a modified Ebers-Moll model. Includes an optional substrate pin and optional
Early voltage effect. 

    # Structural Parameters
        - `use_substrate`: If `true`, a substrate pin connector is available. If `false` it is 
        assumed the substrate is connected to the collector pin.

        - `use_Early`: If `true`, the Early effect is modeled, which takes in to account the effect 
        collector-base voltage variations have on the collector-base depletion region. In many cases this
        effectively means that the collector current has a dependency on the collector-emitter voltage.

        - `use_advanced_continuation`: When false, the `C_jC` and `C_jE` non-linear capacitance curves use 
        a simplified linear continuation starting when `V_BC` and `V_BE` are 0, respectively. If `true`, the `Z_C` and `Z_E` parameters 
        are used to start the linear continuation at `Phi_C - Z_C` and `Phi_E - Z_E`. 

    # Connectors
        - `b` Base Pin
        - `c` Collector Pin
        - `e` Emitter Pin
        - `s` Substrate Pin, only available when `use_substrate = true`

    # Parameters
        - `B_F`: Forward beta
        - `B_R`: Reverse beta
        - `Is`: Saturation current
        - `V_T`: Thermal voltage at 300K
        - `V_A`: Inverse Early voltage
        - `phi_C`: Collector junction exponent
        - `phi_E`: Emitter junction exponent
        - `Z_C`: Collector junction offset
        - `Z_E`: Emitter junction offset 
        - `Tau_f`: Forward transit time
        - `Tau_r`: Reverse transit time
        - `C_jC0`: Collector junction capacitance coefficient
        - `C_jE0`: Emitter junction capacitance coefficient
        - `C_CS`: Collector-substrate capacitance
        - `gamma_C`: Collector junction exponent
        - `gamma_E`: Emitter junction exponent
        - `NF`: Forward emission coefficient
        - `NR`: Reverse emission coefficient
"""
@mtkmodel NPN begin
    @variables begin
        V_BE(t)
        V_BC(t)
        ICC(t)
        IEC(t)

        C_jC(t)
        C_jE(t)
        C_DC(t)
        C_DE(t)

        I_sub(t)
        V_sub(t)
        V_CS(t)
    end

    @structural_parameters begin
        use_substrate = false
        use_Early = true
        use_advanced_continuation = false
    end

    @components begin
        b = Pin()
        e = Pin()
        c = Pin()

        if use_substrate
            s = Pin()
        end
    end

    @parameters begin
        B_F = 50.0, [description = "Forward beta"]
        B_R = 0.1, [description = "Reverse beta"]
        Is = 1e-16, [description = "Saturation current"]
        V_T = 0.026, [description = "Thermal voltage at 300K"]

        if use_Early
            V_A = 0.02, [description = "Inverse Early voltage"]
        end

        phi_C = 0.8, [description = "Collector junction scaling factor"]
        phi_E = 0.6, [description = "Emitter junction scaling factor"]

        if use_advanced_continuation
            Z_C = 0.1, [description = "Collector junction offset"]
            Z_E = 0.1, [description = "Emitter junction offset"]
        end

        Tau_f = 0.12e-9, [description = "Forward transit time"]
        Tau_r = 5e-9, [description = "Reverse transit time"]

        C_jC0 = 0.5e-12, [description = "Collector-junction capacitance coefficient"]
        C_jE0 = 0.4e-12, [description = "Emitter-junction capacitance coefficient"]

        C_CS = 1e-12, [description = "Collector-substrate capacitance"]

        gamma_C = 0.5, [description = "Collector junction exponent"]
        gamma_E = 1.0 / 3.0, [description = "Emitter junction exponent"]

        NF = 1.0, [description = "Forward ideality exponent"]
        NR = 1.0, [description = "Reverse ideality exponent"]
    end

    @equations begin
        V_BE ~ b.v - e.v
        V_BC ~ b.v - c.v

        ICC ~ Is * (exp(V_BE / V_T) - 1)
        IEC ~ Is * (exp(V_BC / V_T) - 1)

        if !use_advanced_continuation
            C_jC ~ ifelse(V_BC / phi_C > 0.0, 1 + gamma_C * V_BC / phi_C,
                (C_jC0) / (1 - V_BC / phi_C)^gamma_C)
            C_jE ~ ifelse(V_BE / phi_E > 0.0, 1 + gamma_E * V_BE / phi_E,
                (C_jE0) / (1 - V_BE / phi_E)^gamma_E)
        end

        if use_advanced_continuation
            C_jC ~ if V_BC > phi_C - Z_C
                ((C_jC0 * gamma_C * (1 - ((phi_C - Z_C) / phi_C))^(-gamma_C - 1)) / phi_C) *
                V_BC -
                ((C_jC0 * gamma_C * (1 - ((phi_C - Z_C) / phi_C))^(-gamma_C - 1)) / phi_C) *
                (phi_C - Z_C) + (C_jC0) / (1 - (phi_C - Z_C) / phi_C)^gamma_C
            else
                (C_jC0) / (1 - V_BC / phi_C)^gamma_C
            end

            C_jE ~ if V_BE > phi_E - Z_E
                ((C_jE0 * gamma_E * (1 - ((phi_E - Z_E) / phi_E))^(-gamma_E - 1)) / phi_E) *
                V_BE -
                ((C_jE0 * gamma_E * (1 - ((phi_E - Z_E) / phi_E))^(-gamma_E - 1)) / phi_E) *
                (phi_E - Z_E) + (C_jE0) / (1 - (phi_E - Z_E) / phi_E)^gamma_E
            else
                (C_jE0) / (1 - V_BE / phi_E)^gamma_E
            end
        end

        C_DE ~ Tau_f * (Is / (NF * V_T)) * exp(V_BE / (NF * V_T))
        C_DC ~ Tau_r * (Is / (NR * V_T)) * exp(V_BC / (NR * V_T))

        if use_substrate
            s.i ~ I_sub
            s.v ~ V_sub
            V_CS ~ c.v - V_sub
        end

        if !use_substrate
            V_sub ~ c.v
        end

        I_sub ~ ifelse(use_substrate, -C_CS * D(V_CS), -C_CS * D(V_sub))

        c.i ~
        (ICC - IEC) * ifelse(use_Early, (1 - V_BC * V_A), 1.0) - IEC / B_R -
        (C_jC + C_DC) * D(V_BC) - I_sub
        b.i ~ IEC / B_R + ICC / B_F + (C_jC + C_DC) * D(V_BC) + (C_jE + C_DE) * D(V_BE)
        e.i ~ -c.i - b.i - I_sub
    end
end

"""
    PNP(;name, B_F, B_R, Is, V_T, V_A, phi_C, phi_E, Z_C, Z_E, Tau_f, Tau_r, C_jC0, C_jE0, C_CS, gamma_C, gamma_E, NF, NR)

Creates a PNP Bipolar Junction Transistor following a modified Ebers-Moll model. Includes an optional substrate pin and optional
Early voltage effect. 

    # Structural Parameters
        - `use_substrate`: If `true`, a substrate pin connector is available. If `false` it is 
        assumed the substrate is connected to the collector pin.

        - `use_Early`: If `true`, the Early effect is modeled, which takes in to account the effect 
        collector-base voltage variations have on the collector-base depletion region. In many cases this
        effectively means that the collector current has a dependency on the collector-emitter voltage.

        - `use_advanced_continuation`: When false, the `C_jC` and `C_jE` non-linear capacitance curves use 
        a simplified linear continuation starting when `V_CB` and `V_EB` are 0, respectively. If `true`, the `Z_C` and `Z_E` parameters 
        are used to start the linear continuation at `Phi_C - Z_C` and `Phi_E - Z_E`. 

    # Connectors
        - `b` Base Pin
        - `c` Collector Pin
        - `e` Emitter Pin
        - `s` Substrate Pin, only available when `use_substrate = true`

    # Parameters
        - `B_F`: Forward beta
        - `B_R`: Reverse beta
        - `Is`: Saturation current
        - `V_T`: Thermal voltage at 300K
        - `V_A`: Inverse Early voltage
        - `phi_C`: Collector junction exponent
        - `phi_E`: Emitter junction exponent
        - `Z_C`: Collector junction offset
        - `Z_E`: Emitter junction offset 
        - `Tau_f`: Forward transit time
        - `Tau_r`: Reverse transit time
        - `C_jC0`: Collector junction capacitance coefficient
        - `C_jE0`: Emitter junction capacitance coefficient
        - `C_CS`: Collector-substrate capacitance
        - `gamma_C`: Collector junction exponent
        - `gamma_E`: Emitter junction exponent
        - `NF`: Forward emission coefficient
        - `NR`: Reverse emission coefficient
"""
@mtkmodel PNP begin
    @variables begin
        V_EB(t)
        V_CB(t)
        ICC(t)
        IEC(t)

        C_jC(t)
        C_jE(t)
        C_DC(t)
        C_DE(t)

        I_sub(t)
        V_sub(t)
        V_CS(t)
    end

    @structural_parameters begin
        use_substrate = false
        use_Early = true
        use_advanced_continuation = false
    end

    @components begin
        b = Pin()
        e = Pin()
        c = Pin()

        if use_substrate
            s = Pin()
        end
    end

    @parameters begin
        B_F = 50.0, [description = "Forward beta"]
        B_R = 0.1, [description = "Reverse beta"]
        Is = 1e-16, [description = "Saturation current"]
        V_T = 0.026, [description = "Thermal voltage at 300K"]

        if use_Early
            V_A = 0.02, [description = "Inverse Early voltage"]
        end

        phi_C = 0.8, [description = "Collector junction scaling factor"]
        phi_E = 0.6, [description = "Emitter junction scaling factor"]

        if use_advanced_continuation
            Z_C = 0.1, [description = "Collector junction offset"]
            Z_E = 0.1, [description = "Emitter junction offset"]
        end

        Tau_f = 0.12e-9, [description = "Forward transit time"]
        Tau_r = 5e-9, [description = "Reverse transit time"]

        C_jC0 = 0.5e-12, [description = "Collector-junction capacitance coefficient"]
        C_jE0 = 0.4e-12, [description = "Emitter-junction capacitance coefficient"]

        C_CS = 1e-12, [description = "Collector-substrate capacitance"]

        gamma_C = 0.5, [description = "Collector junction exponent"]
        gamma_E = 1.0 / 3.0, [description = "Emitter junction exponent"]

        NF = 1.0, [description = "Forward ideality exponent"]
        NR = 1.0, [description = "Reverse ideality exponent"]
    end

    @equations begin
        V_EB ~ e.v - b.v
        V_CB ~ c.v - b.v

        ICC ~ Is * (exp(V_EB / V_T) - 1)
        IEC ~ Is * (exp(V_CB / V_T) - 1)

        if !use_advanced_continuation
            C_jC ~ ifelse(V_CB / phi_C > 0.0, 1 + gamma_C * V_CB / phi_C,
                (C_jC0) / (1 - V_CB / phi_C)^gamma_C)
            C_jE ~ ifelse(V_EB / phi_E > 0.0, 1 + gamma_E * V_EB / phi_E,
                (C_jE0) / (1 - V_EB / phi_E)^gamma_E)
        end

        if use_advanced_continuation
            C_jC ~ if V_CB > phi_C - Z_C
                ((C_jC0 * gamma_C * (1 - ((phi_C - Z_C) / phi_C))^(-gamma_C - 1)) / phi_C) *
                V_CB -
                ((C_jC0 * gamma_C * (1 - ((phi_C - Z_C) / phi_C))^(-gamma_C - 1)) / phi_C) *
                (phi_C - Z_C) + (C_jC0) / (1 - (phi_C - Z_C) / phi_C)^gamma_C
            else
                (C_jC0) / (1 - V_CB / phi_C)^gamma_C
            end

            C_jE ~ if V_EB > phi_E - Z_E
                ((C_jE0 * gamma_E * (1 - ((phi_E - Z_E) / phi_E))^(-gamma_E - 1)) / phi_E) *
                V_EB -
                ((C_jE0 * gamma_E * (1 - ((phi_E - Z_E) / phi_E))^(-gamma_E - 1)) / phi_E) *
                (phi_E - Z_E) + (C_jE0) / (1 - (phi_E - Z_E) / phi_E)^gamma_E
            else
                (C_jE0) / (1 - V_EB / phi_E)^gamma_E
            end
        end

        C_DE ~ Tau_f * (Is / (NF * V_T)) * exp(V_EB / (NF * V_T))
        C_DC ~ Tau_r * (Is / (NR * V_T)) * exp(V_CB / (NR * V_T))

        if use_substrate
            s.i ~ I_sub
            s.v ~ V_sub
            V_CS ~ c.v - V_sub
        end

        if !use_substrate
            V_sub ~ c.v
        end

        I_sub ~ ifelse(use_substrate, -C_CS * D(V_CS), -C_CS * D(V_sub))

        c.i ~
        IEC / B_R - (ICC - IEC) * ifelse(use_Early, (1 - V_CB * V_A), 1.0) +
        (C_jC + C_DC) * D(V_CB) - I_sub
        b.i ~ -IEC / B_R - ICC / B_F - (C_jC + C_DC) * D(V_CB) - (C_jE + C_DE) * D(V_EB)
        e.i ~ -c.i - b.i - I_sub
    end
end
