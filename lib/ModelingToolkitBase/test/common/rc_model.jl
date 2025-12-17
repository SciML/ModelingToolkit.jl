import ModelingToolkitStandardLibrary.Electrical as El
import ModelingToolkitStandardLibrary.Blocks as Bl
using SciCompDSL

@component function RCModel(; name)
    pars = @parameters begin
        R = 1.0
        C = 1.0
        V = 1.0
    end
    systems = @named begin
        resistor = El.Resistor(R = R)
        capacitor = El.Capacitor(C = C)
        shape = Bl.Constant(k = V)
        source = El.Voltage()
        ground = El.Ground()
    end
    eqs = [
        connect(shape.output, source.V)
        connect(source.p, resistor.p)
        connect(resistor.n, capacitor.p)
        connect(capacitor.n, source.n)
        connect(capacitor.n, ground.g)
    ]
    System(eqs, t, [], pars; systems, name)
end

@named rc_model = RCModel()
