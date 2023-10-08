using ModelingToolkit, Test
using ModelingToolkit: get_gui_metadata,
    VariableDescription, getdefault, RegularConnector, get_ps, getname
using URIs: URI
using Distributions
using Unitful

ENV["MTK_ICONS_DIR"] = "$(@__DIR__)/icons"

@testset "Comprehensive Test of Parsing Models (with an RC Circuit)" begin
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

    @named capacitor = Capacitor(C = 10, v = 10.0)
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
        @structural_parameters begin
            R_val = 10
            C_val = 10
            k_val = 10
        end
        @components begin
            resistor = Resistor(; R = R_val)
            capacitor = Capacitor(; C = C_val)
            source = Voltage()
            constant = Constant(; k = k_val)
            ground = Ground()
        end

        @equations begin
            connect(constant.output, source.V)
            connect(source.p, resistor.p)
            connect(resistor.n, capacitor.p)
            connect(capacitor.n, source.n, ground.g)
        end
    end

    C_val = 20
    R_val = 20
    res__R = 100
    @mtkbuild rc = RC(; C_val, R_val, resistor.R = res__R)
    resistor = getproperty(rc, :resistor; namespace = false)
    @test getname(rc.resistor) === getname(resistor)
    @test getname(rc.resistor.R) === getname(resistor.R)
    @test getname(rc.resistor.v) === getname(resistor.v)
    # Test that `resistor.R` overrides `R_val` in the argument.
    @test getdefault(rc.resistor.R) == res__R != R_val
    # Test that `C_val` passed via argument is set as default of C.
    @test getdefault(rc.capacitor.C) == C_val
    # Test that `k`'s default value is unchanged.
    @test getdefault(rc.constant.k) == RC.structure[:kwargs][:k_val]
    @test getdefault(rc.capacitor.v) == 0.0

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

    @test length(equations(rc)) == 1
end

@testset "Parameters and Structural parameters in various modes" begin
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
        @structural_parameters begin
            l = 1
            func
        end
    end

    kval = 5
    @named model = MockModel(; kval, cval = 1, func = identity)

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
end

@testset "Defaults of subcomponents MTKModel" begin
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
    params = get_ps(a)
    @test isequal(getdefault(a.b.i), params[1])
    @test isequal(getdefault(a.b.j), 1 / params[1])
    @test getdefault(a.b.k) == 1

    @named a = A(p = 10, b.i = 20, b.j = 30, b.k = 40)
    @test getdefault(a.b.i) == 20
    @test getdefault(a.b.j) == 30
    @test getdefault(a.b.k) == 40
end

@testset "Metadata in variables" begin
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
end

@testset "Connector with parameters, equations..." begin
    @connector A begin
        @extend (e,) = extended_e = E()
        @icon "pin.png"
        @parameters begin
            p
        end
        @variables begin
            v(t)
        end
        @components begin
            cc = C()
        end
        @equations begin
            e ~ 0
        end
    end

    @connector C begin
        c(t)
    end

    @connector E begin
        e(t)
    end

    @named aa = A()
    @test aa.connector_type == RegularConnector()

    @test A.isconnector == true

    @test A.structure[:parameters] == Dict(:p => Dict())
    @test A.structure[:extend] == [[:e], :extended_e, :E]
    @test A.structure[:equations] == ["e ~ 0"]
    @test A.structure[:kwargs] == Dict(:p => nothing, :v => nothing)
    @test A.structure[:components] == [[:cc, :C]]
end
