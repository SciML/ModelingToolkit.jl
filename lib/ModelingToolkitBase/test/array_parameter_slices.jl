using ModelingToolkitBase, Symbolics, SymbolicIndexingInterface, Test

@testset "Parameter defaults and bindings preserve array slices" begin
    @independent_variables time
    @parameters parent[1:2, 1:3]
    values = [1.0 3.0 5.0; 2.0 4.0 6.0]
    for indices in ((Colon(), 1:1), (1:2, 1:2:3), (2:-1:1, 3:-2:1))
        expected = values[indices...]
        n, m = size(expected)
        @parameters selected[1:n, 1:m]
        for bound in (false, true)
            defaults = Pair[parent => values]
            bindings = Pair[]
            push!(bound ? bindings : defaults, selected => parent[indices...])
            model = complete(
                System(
                    Equation[], time, [], [parent, selected];
                    name = :parameters_only, initial_conditions = defaults, bindings
                )
            )
            parameters = MTKParameters(model, [])
            @test getp(model, selected)(parameters) == expected
            overridden = MTKParameters(model, [parent => 2values])
            @test getp(model, selected)(overridden) == 2expected
        end
    end
end
