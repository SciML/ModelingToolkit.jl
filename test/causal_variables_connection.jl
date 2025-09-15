using ModelingToolkit, ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit: t_nounits as t, D_nounits as D

@testset "Error checking" begin
    @variables begin
        x(t)
        y(t), [input = true]
        z(t), [output = true]
        w(t)
        v(t), [input = true]
        u(t), [output = true]
        xarr(t)[1:4], [output = true]
        yarr(t)[1:2, 1:2], [input = true]
    end
    @parameters begin
        p, [input = true]
        q, [output = true]
    end

    @test_throws ["p", "kind", "VARIABLE", "PARAMETER"] connect(z, p)
    @test_throws ["q", "kind", "VARIABLE", "PARAMETER"] connect(q, y)
    @test_throws ["p", "kind", "VARIABLE", "PARAMETER"] connect(z, y, p)

    @test_throws ["unique"] connect(z, y, y)

    @test_throws ["same size"] connect(xarr, yarr)

    @test_throws ArgumentError connect(x, y)
end

@testset "Connection expansion" begin
    @named P = FirstOrder(k = 1, T = 1)
    @named C = Gain(; k = -1)

    eqs = [connect(P.output.u, C.input.u)
           connect(C.output.u, P.input.u)]
    sys1 = System(eqs, t, systems = [P, C], name = :hej)
    sys = expand_connections(sys1)
    @test any(isequal(C.input.u ~ P.output.u), equations(sys))
    @test any(isequal(P.input.u ~ C.output.u), equations(sys))

    @named sysouter = System(Equation[], t; systems = [sys1])
    sys = expand_connections(sysouter)
    @test any(isequal(sys1.C.input.u ~ sys1.P.output.u), equations(sys))
    @test any(isequal(sys1.P.input.u ~ sys1.C.output.u), equations(sys))
end

@testset "With Analysis Points" begin
    @named P = FirstOrder(k = 1, T = 1)
    @named C = Gain(; k = -1)

    ap = AnalysisPoint(:plant_input)
    eqs = [connect(P.output, C.input), connect(C.output.u, ap, P.input.u)]
    sys = System(eqs, t, systems = [P, C], name = :hej)
    @named nested_sys = System(Equation[], t; systems = [sys])

    test_cases = [
        ("inner", sys, sys.plant_input),
        ("nested", nested_sys, nested_sys.hej.plant_input),
        ("inner - Symbol", sys, :plant_input),
        ("nested - Symbol", nested_sys, nameof(sys.plant_input))
    ]

    @testset "get_sensitivity - $name" for (name, sys, ap) in test_cases
        matrices, _ = get_sensitivity(sys, ap)
        @test matrices.A[] == -2
        @test matrices.B[] * matrices.C[] == -1 # either one negative
        @test matrices.D[] == 1
    end

    @testset "get_comp_sensitivity - $name" for (name, sys, ap) in test_cases
        matrices, _ = get_comp_sensitivity(sys, ap)
        @test matrices.A[] == -2
        @test matrices.B[] * matrices.C[] == 1 # both positive or negative
        @test matrices.D[] == 0
    end

    @testset "get_looptransfer - $name" for (name, sys, ap) in test_cases
        matrices, _ = get_looptransfer(sys, ap)
        @test matrices.A[] == -1
        @test matrices.B[] * matrices.C[] == -1 # either one negative
        @test matrices.D[] == 0
    end

    @testset "open_loop - $name" for (name, sys, ap) in test_cases
        open_sys, (du, u) = open_loop(sys, ap)
        matrices, _ = linearize(open_sys, [du], [u])
        @test matrices.A[] == -1
        @test matrices.B[] * matrices.C[] == -1 # either one negative
        @test matrices.D[] == 0
    end
end

@testset "Outside input to inside input connection" begin
    @mtkmodel Inner begin
        @variables begin
            x(t), [input = true]
            y(t), [output = true]
        end
        @equations begin
            y ~ x
        end
    end
    @mtkmodel Outer begin
        @variables begin
            u(t), [input = true]
            v(t), [output = true]
        end
        @components begin
            inner = Inner()
        end
        @equations begin
            connect(u, inner.x)
            connect(inner.y, v)
        end
    end
    @named sys = Outer()
    ss = toggle_namespacing(sys, false)
    eqs = equations(expand_connections(sys))
    @test issetequal(eqs, [ss.inner.x ~ ss.u
                           ss.inner.y ~ ss.inner.x
                           ss.v ~ ss.inner.y])
end
