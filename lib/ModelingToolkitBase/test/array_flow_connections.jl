using ModelingToolkitBase, Symbolics, Test
using ModelingToolkitBase: t_nounits as t

@connector function ArrayFlowPort(; name, rows = 1, columns = 1)
    vars = @variables potential(t)[1:rows, 1:columns] flow(t)[1:rows, 1:columns] [connect = Flow]
    return System(Equation[], t, vars, []; name)
end

function connection_block_size(x)
    x = Symbolics.value(x)
    if Symbolics.iscall(x)
        return 1 + sum(connection_block_size, Symbolics.arguments(x); init = 0)
    elseif x isa AbstractArray
        return 1 + sum(connection_block_size, x; init = 0)
    end
    return 1
end

@testset "Array flow connections retain equation blocks" begin
    counts = Int[]
    sizes = Int[]
    for rows in (1, 3, 100), columns in (1, 2)
        @named a = ArrayFlowPort(; rows, columns)
        @named b = ArrayFlowPort(; rows, columns)
        source = System(
            [a.potential ~ fill(300.0, rows, columns)], t;
            systems = [a], name = :source
        )
        sink = System(
            [b.flow ~ 2 .* (b.potential .- 298.15)], t;
            systems = [b], name = :sink
        )
        network = System(
            [connect(source.a, sink.b)], t;
            systems = [source, sink], name = :network
        )
        sys, info = expand_connections(network, Val(true))
        eqs = equations(sys)
        ceqs = eqs[info.is_connection_equation]
        push!(counts, length(ceqs))
        push!(sizes, sum(eq -> connection_block_size(eq.lhs) + connection_block_size(eq.rhs), ceqs))
        point = Dict{Any, Any}()
        for (field, value) in (
                (source.a.potential, 300.0), (sink.b.potential, 300.0),
                (source.a.flow, -3.7), (sink.b.flow, 3.7),
            )
            for v in Symbolics.scalarize(field)
                point[Symbolics.value(v)] = value
            end
        end
        for eq in eqs
            scalar_eqs = Symbolics.scalarize(eq)
            for e in (scalar_eqs isa AbstractArray ? scalar_eqs : (scalar_eqs,))
                @test abs(Symbolics.value(Symbolics.substitute(e.lhs - e.rhs, point))) < 1.0e-12
            end
        end
    end
    @test allequal(counts)
    @test allequal(sizes)
end
