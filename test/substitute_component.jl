using ModelingToolkit, ModelingToolkitStandardLibrary, Test
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
using OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D, renamespace,
    NAMESPACE_SEPARATOR as NS

@component function SignalInterface(; name)
    pars = @parameters begin
    end

    systems = @named begin
        output = RealOutput()
    end

    vars = @variables begin
    end

    equations = Equation[]

    return System(equations, t, vars, pars; name, systems)
end

@component function TwoComponent(; name)
    pars = @parameters begin
    end

    systems = @named begin
        component1 = OnePort()
        component2 = OnePort()
        source = Voltage()
        signal = SignalInterface()
        ground = Ground()
    end

    vars = @variables begin
    end

    equations = Equation[
        connect(signal.output.u, source.V.u),
        connect(source.p, component1.p),
        connect(component1.n, component2.p),
        connect(component2.n, source.n, ground.g),
    ]

    return System(equations, t, vars, pars; name, systems)
end

@component function RC(; name, R = 1.0, C = 1.0, V = 1.0)
    pars = @parameters begin
        R = R
        C = C
        V = V
    end

    systems = @named begin
        component1 = Resistor(R = R)
        component2 = Capacitor(C = C, v = 0.0)
        source = Voltage()
        constant = Constant(k = V)
        ground = Ground()
    end

    vars = @variables begin
    end

    equations = Equation[
        connect(constant.output, source.V),
        connect(source.p, component1.p),
        connect(component1.n, component2.p),
        connect(component2.n, source.n, ground.g),
    ]

    return System(equations, t, vars, pars; name, systems)
end

@testset "Replacement with connections works" begin
    @named templated = TwoComponent()
    @named component1 = Resistor(R = 1.0)
    @named component2 = Capacitor(C = 1.0, v = 0.0)
    @named signal = Constant(k = 1.0)
    rsys = substitute_component(templated, templated.component1 => component1)
    rcsys = substitute_component(rsys, rsys.component2 => component2)
    rcsys = substitute_component(rcsys, rcsys.signal => signal)

    @named reference = RC()

    sys1 = mtkcompile(rcsys)
    sys2 = mtkcompile(reference)
    @test isequal(unknowns(sys1), unknowns(sys2))
    @test isequal(equations(sys1), equations(sys2))

    prob1 = ODEProblem(sys1, [], (0.0, 10.0))
    prob2 = ODEProblem(sys2, [], (0.0, 10.0))

    sol1 = solve(prob1, Tsit5())
    sol2 = solve(prob2, Tsit5(); saveat = sol1.t)
    @test sol1.u ≈ sol2.u atol = 1.0e-8
end

@component function BadOnePort1(; name)
    pars = @parameters begin
    end

    systems = @named begin
        p = Pin()
        n = Pin()
    end

    vars = @variables begin
        i(t)
    end

    equations = Equation[
        0 ~ p.i + n.i,
        i ~ p.i,
    ]

    return System(equations, t, vars, pars; name, systems)
end

@connector function BadPin1(; name, v = nothing)
    vars = @variables begin
        v(t) = v
    end
    System(Equation[], t, vars, []; name)
end

@component function BadOnePort2(; name)
    pars = @parameters begin
    end

    systems = @named begin
        p = BadPin1()
        n = Pin()
    end

    vars = @variables begin
        v(t)
        i(t)
    end

    equations = Equation[
        0 ~ p.v + n.v,
        v ~ p.v,
    ]

    return System(equations, t, vars, pars; name, systems)
end

@connector function BadPin2(; name, v = nothing, i = nothing)
    vars = @variables begin
        v(t) = v
        i(t) = i
    end
    System(Equation[], t, vars, []; name)
end

@component function BadOnePort3(; name)
    pars = @parameters begin
    end

    systems = @named begin
        p = BadPin2()
        n = Pin()
    end

    vars = @variables begin
        v(t)
        i(t)
    end

    equations = Equation[
        0 ~ p.v + n.v,
        v ~ p.v,
    ]

    return System(equations, t, vars, pars; name, systems)
end

@connector function BadPin3(; name, v = nothing, i = nothing)
    vars = @variables begin
        v(t) = v, [input = true]
        i(t) = i, [connect = Flow]
    end
    System(Equation[], t, vars, []; name)
end

@component function BadOnePort4(; name)
    pars = @parameters begin
    end

    systems = @named begin
        p = BadPin3()
        n = Pin()
    end

    vars = @variables begin
        v(t)
        i(t)
    end

    equations = Equation[
        0 ~ p.v + n.v,
        v ~ p.v,
    ]

    return System(equations, t, vars, pars; name, systems)
end

@connector function BadPin4(; name, v = nothing, i = nothing)
    vars = @variables begin
        v(t) = v, [output = true]
        i(t) = i, [connect = Flow]
    end
    System(Equation[], t, vars, []; name)
end

@component function BadOnePort5(; name)
    pars = @parameters begin
    end

    systems = @named begin
        p = BadPin4()
        n = Pin()
    end

    vars = @variables begin
        v(t)
        i(t)
    end

    equations = Equation[
        0 ~ p.v + n.v,
        v ~ p.v,
    ]

    return System(equations, t, vars, pars; name, systems)
end

@component function BadPin5(; name)
    pars = @parameters begin
    end

    systems = @named begin
    end

    vars = @variables begin
        v(t)
        i(t), [connect = Flow]
    end

    equations = Equation[]

    return System(equations, t, vars, pars; name, systems)
end

@component function BadOnePort6(; name)
    pars = @parameters begin
    end

    systems = @named begin
        p = BadPin5()
        n = Pin()
    end

    vars = @variables begin
        v(t)
        i(t)
    end

    equations = Equation[
        0 ~ p.v + n.v,
        v ~ p.v,
    ]

    return System(equations, t, vars, pars; name, systems)
end

@connector function BadPin6(; name, i = nothing)
    vars = @variables begin
        i(t) = i, [connect = Flow]
    end
    System(Equation[], t, vars, []; name)
end

@component function BadOnePort7(; name)
    pars = @parameters begin
    end

    systems = @named begin
        p = BadPin6()
        n = Pin()
    end

    vars = @variables begin
        v(t)
        i(t)
    end

    equations = Equation[
        0 ~ p.i + n.i,
        i ~ p.i,
    ]

    return System(equations, t, vars, pars; name, systems)
end

@component function BadOnePort8(; name)
    pars = @parameters begin
    end

    systems = @named begin
        n = Pin()
    end

    vars = @variables begin
        v(t)
        i(t)
    end

    equations = Equation[]

    return System(equations, t, vars, pars; name, systems)
end

@testset "Error checking" begin
    @named templated = TwoComponent()
    @named component1 = Resistor(R = 1.0)
    @named component2 = Capacitor(C = 1.0, v = 0.0)
    @test_throws ["LHS", "cannot be completed"] substitute_component(
        templated, complete(templated.component1) => component1
    )
    @test_throws ["RHS", "cannot be completed"] substitute_component(
        templated, templated.component1 => complete(component1)
    )
    @test_throws ["RHS", "not be namespaced"] substitute_component(
        templated, templated.component1 => renamespace(templated, component1)
    )
    @named resistor = Resistor(R = 1.0)
    @test_throws ["RHS", "same name"] substitute_component(
        templated, templated.component1 => resistor
    )

    @testset "Different indepvar" begin
        @independent_variables tt
        @named empty = System(Equation[], t)
        @named outer = System(Equation[], t; systems = [empty])
        @named empty = System(Equation[], tt)
        @test_throws ["independent variable"] substitute_component(
            outer, outer.empty => empty
        )
    end

    @named component1 = BadOnePort1()
    @test_throws ["RHS", "unknown", "v(t)"] substitute_component(
        templated, templated.component1 => component1
    )

    @named component1 = BadOnePort2()
    @test_throws ["component1$(NS)p", "i(t)"] substitute_component(
        templated, templated.component1 => component1
    )

    @named component1 = BadOnePort3()
    @test_throws ["component1$(NS)p$(NS)i", "Flow"] substitute_component(
        templated, templated.component1 => component1
    )

    @named component1 = BadOnePort4()
    @test_throws ["component1$(NS)p$(NS)v", "differing causality", "input"] substitute_component(
        templated, templated.component1 => component1
    )

    @named component1 = BadOnePort5()
    @test_throws ["component1$(NS)p$(NS)v", "differing causality", "output"] substitute_component(
        templated, templated.component1 => component1
    )

    @named component1 = BadOnePort6()
    @test_throws ["templated$(NS)component1$(NS)p", "not a connector"] substitute_component(
        templated, templated.component1 => component1
    )

    @named component1 = BadOnePort7()
    @test_throws ["templated$(NS)component1$(NS)p", "DomainConnector", "RegularConnector"] substitute_component(
        templated, templated.component1 => component1
    )

    @named component1 = BadOnePort8()
    @test_throws ["templated$(NS)component1", "subsystem p"] substitute_component(
        templated, templated.component1 => component1
    )
end
