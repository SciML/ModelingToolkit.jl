using ModelingToolkit, Test
using ModelingToolkit: get_gui_metadata
using URIs: URI

ENV["MTK_ICONS_DIR"] = "$(@__DIR__)/icons"

@connector RealInput begin
    u(t), [input = true]
end
@connector RealOutput begin
    u(t), [output = true]
end
@model Constant begin
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

@variables t
D = Differential(t)

@connector Pin begin
    v(t) = 0                  # Potential at the pin [V]
    i(t), [connect = Flow]    # Current flowing into the pin [A]
    @icon "pin.png"
end

@model OnePort begin
    @components begin
        p = Pin()
        n = Pin()
    end
    @variables begin
        v(t)
        i(t)
    end
    @icon "oneport.png"
    @equations begin
        v ~ p.v - n.v
        0 ~ p.i + n.i
        i ~ p.i
    end
end

@model Ground begin
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
@model Resistor begin
    @extend v, i = oneport = OnePort()
    @parameters begin
        R = 1
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

@model Capacitor begin
    @extend v, i = oneport = OnePort()
    @parameters begin
        C = 1
    end
    @icon "https://upload.wikimedia.org/wikipedia/commons/7/78/Capacitor_symbol.svg"
    @equations begin
        D(v) ~ i / C
    end
end

@model Voltage begin
    @extend v, i = oneport = OnePort()
    @components begin
        V = RealInput()
    end
    @equations begin
        v ~ V.u
    end
end

@model RC begin
    @components begin
        resistor = Resistor()
        capacitor = Capacitor()
        source = Voltage()
        constant = Constant()
        ground = Ground()
    end
    @equations begin
        connect(constant.output, source.V)
        connect(source.p, resistor.p)
        connect(resistor.n, capacitor.p)
        connect(capacitor.n, source.n, ground.g)
    end
end
@named rc = RC()

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
