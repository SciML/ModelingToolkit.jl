"""
    Force(; name, use_support = false)

Input signal acting as external force on a flange
"""
@mtkmodel Force begin
    @extend (flange,) = partial_element = PartialElementaryOneFlangeAndSupport2(;
        use_support = false)
    @parameters begin
        s = 0
    end
    @components begin
        f = RealInput()
    end
    @equations begin
        flange.f ~ -f.u
    end
end
