import ModelingToolkitStandardLibrary.Electrical as El
import ModelingToolkitStandardLibrary.Blocks as Bl

@mtkmodel RCModel begin
    @parameters begin
        R = 1.0
        C = 1.0
        V = 1.0
    end
    @components begin
        resistor = El.Resistor(R = R)
        capacitor = El.Capacitor(C = C)
        shape = Bl.Constant(k = V)
        source = El.Voltage()
        ground = El.Ground()
    end
    @equations begin
        connect(shape.output, source.V)
        connect(source.p, resistor.p)
        connect(resistor.n, capacitor.p)
        connect(capacitor.n, source.n)
        connect(capacitor.n, ground.g)
    end
end

@named rc_model = RCModel()
