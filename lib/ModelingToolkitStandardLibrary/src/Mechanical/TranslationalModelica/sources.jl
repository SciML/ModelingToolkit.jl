"""
    Force(; name, use_support = false)

Input signal acting as external force on a flange
"""
@mtkmodel Force begin
    @extend (flange,) = partial_element = PartialElementaryOneFlangeAndSupport2(;
        use_support = false)
    @components begin
        f = RealInput() # Accelerating force acting at flange (= -flange.tau)
    end
    @equations begin
        flange.f ~ -f.u
    end
end

"""
    Position(; name, exact = false, f_crit = 50)

Forced movement of a flange according to a reference position

The input signal `s_ref` defines the reference position in [m]. Flange flange is forced to move relative to the support connector according to this reference motion. According to parameter `exact`, this is done in the following way:

- `exact=true`: The reference position is treated exactly. This is only possible, if the input signal is defined by an analytical function which can be differentiated at least twice. If this prerequisite is fulfilled, the Modelica translator will differentiate the input signal twice in order to compute the reference acceleration of the flange.
- `exact=false`: The reference position is filtered and the second derivative of the filtered curve is used to compute the reference acceleration of the flange. This second derivative is not computed by numerical differentiation but by an appropriate realization of the filter. For filtering, a second order Bessel filter is used. The critical frequency (also called cut-off frequency) of the filter is defined via parameter `f_crit` in [Hz]. This value should be selected in such a way that it is higher as the essential low frequencies in the signal.

The input signal can be provided from one of the signal generator blocks of the block library `Blocks.Sources`.
"""
@mtkmodel Position begin
    @extend (s,) = ptf = PartialElementaryOneFlangeAndSupport2()
    @structural_parameters begin
        exact = false
    end
    @parameters begin
        f_crit = 50
    end
    @variables begin
        v(t)
        a(t)
    end
    @components begin
        s_ref = RealInput()
    end
    begin
        w_crit = 2π * f_crit
        af = 1.3617
        bf = 0.6180
    end
    @equations begin
        if exact
            s ~ s_ref.u
        else
            a ~ ((s_ref.u - s) * w_crit - af * v) * (w_crit / bf)
        end
        v ~ D(s)
        a ~ D(v)
    end
end
