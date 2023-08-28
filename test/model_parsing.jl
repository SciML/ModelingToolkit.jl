using ModelingToolkit, Test
using ModelingToolkit: get_gui_metadata, VariableDescription, getdefault
using URIs: URI
using Distributions
using Unitful

ENV["MTK_ICONS_DIR"] = "$(@__DIR__)/icons"

@connector RealInput begin
    u(t), [input = true, unit = u"V"]
end
@connector RealOutput begin
    u(t), [output = true, unit = u"V"]
end
@mtkmodel Constant begin
    @components begin
        output = RealOutput()
    end
    @parameters begin
        k, [description = "Constant output value of block"]
    end
    @equations begin
        output.u ~ k
    end
end

@variables t [unit = u"s"]
D = Differential(t)

@connector Pin begin
    v(t), [unit = u"V"]                    # Potential at the pin [V]
    i(t), [connect = Flow, unit = u"A"]    # Current flowing into the pin [A]
    @icon "pin.png"
end

@named p = Pin(; v = π)
@test getdefault(p.v) == π
@test Pin.isconnector == true

@mtkmodel OnePort begin
    @components begin
        p = Pin()
        n = Pin()
    end
    @variables begin
        v(t), [unit = u"V"]
        i(t), [unit = u"A"]
    end
    @icon "oneport.png"
    @equations begin
        v ~ p.v - n.v
        0 ~ p.i + n.i
        i ~ p.i
    end
end

@test OnePort.isconnector == false

@mtkmodel Ground begin
    @components begin
        g = Pin()
    end
    @icon begin
        read(abspath(ENV["MTK_ICONS_DIR"], "ground.svg"), String)
    end
    @equations begin
        g.v ~ 0
    end
end

resistor_log = "$(@__DIR__)/logo/resistor.svg"
@mtkmodel Resistor begin
    @extend v, i = oneport = OnePort()
    @parameters begin
        R, [unit = u"Ω"]
    end
    @icon begin
        """<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="80" height="30">
  <path d="M10 15
l15 0
l2.5 -5
l5 10
l5 -10
l5 10
l5 -10
l5 10
l2.5 -5
l15 0" stroke="black" stroke-width="1" stroke-linejoin="bevel" fill="none"></path>
</svg>
"""
    end
    @equations begin
        v ~ i * R
    end
end

@mtkmodel Capacitor begin
    @parameters begin
        C, [unit = u"F"]
    end
    @extend v, i = oneport = OnePort(; v = 0.0)
    @icon "https://upload.wikimedia.org/wikipedia/commons/7/78/Capacitor_symbol.svg"
    @equations begin
        D(v) ~ i / C
    end
end

@named capacitor = Capacitor(C = 10, oneport.v = 10.0)
@test getdefault(capacitor.v) == 10.0

@mtkmodel Voltage begin
    @extend v, i = oneport = OnePort()
    @components begin
        V = RealInput()
    end
    @equations begin
        v ~ V.u
    end
end

@mtkmodel RC begin
    @components begin
        resistor = Resistor(; R)
        capacitor = Capacitor(; C = 10)
        source = Voltage()
        constant = Constant(; k = 1)
        ground = Ground()
    end
    @equations begin
        connect(constant.output, source.V)
        connect(source.p, resistor.p)
        connect(resistor.n, capacitor.p)
        connect(capacitor.n, source.n, ground.g)
    end
end

@named rc = RC(; resistor.R = 20)
@test getdefault(rc.resistor.R) == 20
@test getdefault(rc.capacitor.C) == 10
@test getdefault(rc.capacitor.v) == 0.0
@test getdefault(rc.constant.k) == 1

@test get_gui_metadata(rc.resistor).layout == Resistor.structure[:icon] ==
      read(joinpath(ENV["MTK_ICONS_DIR"], "resistor.svg"), String)
@test get_gui_metadata(rc.ground).layout ==
      read(abspath(ENV["MTK_ICONS_DIR"], "ground.svg"), String)
@test get_gui_metadata(rc.capacitor).layout ==
      URI("https://upload.wikimedia.org/wikipedia/commons/7/78/Capacitor_symbol.svg")
@test OnePort.structure[:icon] ==
      URI("file:///" * abspath(ENV["MTK_ICONS_DIR"], "oneport.png"))
@test ModelingToolkit.get_gui_metadata(rc.resistor.p).layout == Pin.structure[:icon] ==
      URI("file:///" * abspath(ENV["MTK_ICONS_DIR"], "pin.png"))

@test length(equations(structural_simplify(rc))) == 1

@mtkmodel MockModel begin
    @parameters begin
        a
        b(t)
        cval
        jval
        kval
        c(t) = cval + jval
        d = 2
        e, [description = "e"]
        f = 3, [description = "f"]
        h(t), [description = "h(t)"]
        i(t) = 4, [description = "i(t)"]
        j(t) = jval, [description = "j(t)"]
        k = kval, [description = "k"]
    end
end

kval = 5
@named model = MockModel(; kval, cval = 1)

@test hasmetadata(model.e, VariableDescription)
@test hasmetadata(model.f, VariableDescription)
@test hasmetadata(model.h, VariableDescription)
@test hasmetadata(model.i, VariableDescription)
@test hasmetadata(model.j, VariableDescription)
@test hasmetadata(model.k, VariableDescription)

model = complete(model)
@test getdefault(model.cval) == 1
@test isequal(getdefault(model.c), model.cval + model.jval)
@test getdefault(model.d) == 2
@test_throws KeyError getdefault(model.e)
@test getdefault(model.f) == 3
@test getdefault(model.i) == 4
@test isequal(getdefault(model.j), model.jval)
@test isequal(getdefault(model.k), model.kval)

@mtkmodel A begin
    @parameters begin
        p
    end
    @components begin
        b = B(i = p, j = 1 / p, k = 1)
    end
end

@mtkmodel B begin
    @parameters begin
        i
        j
        k
    end
end

@named a = A(p = 10)
getdefault(a.b.i) == 10
getdefault(a.b.j) == 0.1
getdefault(a.b.k) == 1

@named a = A(p = 10, b.i = 20, b.j = 30, b.k = 40)
getdefault(a.b.i) == 20
getdefault(a.b.j) == 30
getdefault(a.b.k) == 40

metadata = Dict(:description => "Variable to test metadata in the Model.structure",
    :input => true, :bounds => (-1, 1), :connection_type => :Flow, :integer => true,
    :binary => false, :tunable => false, :disturbance => true, :dist => Normal(1, 1))

@connector MockMeta begin
    m(t),
    [description = "Variable to test metadata in the Model.structure",
        input = true, bounds = (-1, 1), connect = Flow, integer = true,
        binary = false, tunable = false, disturbance = true, dist = Normal(1, 1)]
end

for (k, v) in metadata
    @test MockMeta.structure[:variables][:m][k] == v
end
