using ModelingToolkit, Test

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
    @equations begin
        g.v ~ 0
    end
end

@model Resistor begin
    @extend v, i = oneport = OnePort()
    @parameters begin
        R = 1
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
@test length(equations(structural_simplify(rc))) == 1
