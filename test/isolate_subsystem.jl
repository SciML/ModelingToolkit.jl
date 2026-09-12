using ModelingToolkit, ModelingToolkitStandardLibrary, Test
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkit: t_nounits as t, get_eqs, get_systems

@testset "Plant extraction between analysis points" begin
    @named plant = FirstOrder(k = 1, T = 1)
    @named controller = Gain(k = -1)
    eqs = [
        connect(controller.output, :plant_input, plant.input)
        connect(plant.output, :plant_output, controller.input)
    ]
    @named closed_loop = System(eqs, t, systems = [plant, controller])
    isolated, input_vars, output_vars = isolate_subsystem(
        closed_loop, :plant_input, :plant_output
    )
    @test isequal(only(input_vars), plant.input.u)
    @test isequal(only(output_vars), plant.output.u)
    @test nameof.(get_systems(isolated)) == [:plant]
end

@testset "multiconnect networks are restricted to the isolated region" begin
    @named plant = FirstOrder(k = 1, T = 1)
    @named inner = Gain(k = 2)
    @named controller = Gain(k = -1)
    @named r1 = Resistor(R = 1.0)
    @named r2 = Resistor(R = 2.0)
    @named gnd = Ground()
    # one network with two disconnected components: the causal plant -> inner edge
    # inside the boundary, and an electrical loop unrelated to it
    network = multiconnect(
        [plant, inner, r1, r2, gnd],
        [
            ConnectionEdge(1, 2, :output, :input),
            ConnectionEdge(3, 4, :n, :p),
            ConnectionEdge(4, 5, :n, :g),
            ConnectionEdge(5, 3, :g, :p),
        ]
    )
    eqs = [
        connect(controller.output, :plant_input, plant.input)
        network
        connect(inner.output, :plant_output, controller.input)
    ]
    @named closed_loop = System(
        eqs, t, systems = [plant, inner, controller, r1, r2, gnd]
    )
    isolated, input_vars, output_vars = isolate_subsystem(
        closed_loop, :plant_input, :plant_output
    )
    @test isequal(only(input_vars), plant.input.u)
    @test isequal(only(output_vars), inner.output.u)
    @test Set(nameof.(get_systems(isolated))) == Set([:plant, :inner])

    conneq = only(filter(eq -> ModelingToolkit.value(eq.lhs) isa Connection, get_eqs(isolated)))
    net = ModelingToolkit.value(conneq.rhs).systems
    @test net isa ConnectionNetwork
    @test nameof.(net.nodes) == [:plant, :inner]
    @test net.edges == [ConnectionEdge(1, 2, :output, :input)]

    expanded = equations(expand_connections(isolated))
    @test any(expanded) do eq
        isequal(eq, inner.input.u ~ plant.output.u) ||
            isequal(eq, plant.output.u ~ inner.input.u)
    end
end
