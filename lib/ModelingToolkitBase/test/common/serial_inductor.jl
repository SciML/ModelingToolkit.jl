import ModelingToolkitStandardLibrary.Electrical as El
import ModelingToolkitStandardLibrary.Blocks as Bl

@mtkmodel LLModel begin
    @components begin
        shape = Bl.Constant(k = 10.0)
        source = El.Voltage()
        resistor = El.Resistor(R = 1.0)
        inductor1 = El.Inductor(L = 1.0e-2)
        inductor2 = El.Inductor(L = 2.0e-2)
        ground = El.Ground()
    end
    @equations begin
        connect(shape.output, source.V)
        connect(source.p, resistor.p)
        connect(resistor.n, inductor1.p)
        connect(inductor1.n, inductor2.p)
        connect(source.n, inductor2.n)
        connect(inductor2.n, ground.g)
    end
end

@named ll_model = LLModel()

@mtkmodel LL2Model begin
    @components begin
        shape = Bl.Constant(k = 10.0)
        source = El.Voltage()
        resistor1 = El.Resistor(R = 1.0)
        resistor2 = El.Resistor(R = 1.0)
        inductor1 = El.Inductor(L = 1.0e-2)
        inductor2 = El.Inductor(L = 2.0e-2)
        ground = El.Ground()
    end
    @equations begin
        connect(shape.output, source.V)
        connect(source.p, inductor1.p)
        connect(inductor1.n, resistor1.p)
        connect(inductor1.n, resistor2.p)
        connect(resistor1.n, resistor2.n)
        connect(resistor2.n, inductor2.p)
        connect(source.n, inductor2.n)
        connect(inductor2.n, ground.g)
    end
end

@named ll2_model = LL2Model()
