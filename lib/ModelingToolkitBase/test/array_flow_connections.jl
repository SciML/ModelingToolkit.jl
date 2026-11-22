using ModelingToolkitBase, Symbolics, SciMLBase, Test
using FillArrays: Zeros
using LinearAlgebra: rank
using ModelingToolkitBase: t_nounits as t, D_nounits as D

@connector function ArrayFlowPort(; name, rows = 1, columns = 1, dynamic = false)
    vars = @variables potential(t)[1:rows, 1:columns] flow(t)[1:rows, 1:columns] [connect = Flow]
    eqs = dynamic ? [D(potential) ~ (flow .- potential)] : Equation[]
    return System(eqs, t, vars, []; name)
end

function connection_block_size(x)
    x = Symbolics.value(x)
    if Symbolics.iscall(x)
        return 1 + sum(connection_block_size, Symbolics.arguments(x); init = 0)
    elseif x isa Zeros
        return 1 + ndims(x)
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

@testset "Signed array flow conservation" begin
    @named a = ArrayFlowPort(; rows = 2, columns = 2)
    @named b = ArrayFlowPort(; rows = 2, columns = 2)
    @named c = ArrayFlowPort(; rows = 2, columns = 2)
    @named outer = ArrayFlowPort(; rows = 2, columns = 2)
    @named inner = System(Equation[], t; systems = [a, b, c])
    @named all_inner = System([connect(inner.a, inner.b, inner.c)], t; systems = [inner])
    @named mixed = System([connect(inner.a, inner.b, outer)], t; systems = [inner, outer])

    a_flow = collect(inner.a.flow)
    b_flow = collect(inner.b.flow)
    c_flow = collect(inner.c.flow)
    outer_flow = collect(outer.flow)
    for (network, expected) in (
            (all_inner, a_flow .+ b_flow .+ c_flow),
            (mixed, a_flow .+ b_flow .- outer_flow),
        )
        eqs = equations(expand_connections(network))
        scalar_eqs = reduce(vcat, Symbolics.scalarize.(eqs))
        for term in expected
            @test (0 ~ term) in scalar_eqs
        end
    end
end

@testset "Array flow port edge cases" begin
    @testset "Open port in a mass-matrix ODE" begin
        @named a = ArrayFlowPort(; rows = 2, columns = 3, dynamic = true)
        @named network = System(Equation[], t; systems = [a])
        p0 = [1.0 3.0 5.0; 2.0 4.0 6.0]
        prob = ODEProblem(
            complete(network), [a.potential => p0, a.flow => zeros(2, 3)],
            (0.0, 1.0); build_initializeprob = false
        )
        @test rank(prob.f.mass_matrix) == 6
        @test prob.f(prob.u0, prob.p, 0.0) == vcat(-vec(p0), zeros(6))
    end

    @testset "Zero-width port in structural compilation" begin
        @variables x(t)
        @named a = ArrayFlowPort(; rows = 0, columns = 3)
        @named network = System([D(x) ~ -x], t, [x], []; systems = [a])
        prob = ODEProblem(mtkcompile(network), [x => 3.0], (0.0, 1.0))
        @test prob.f(prob.u0, prob.p, 0.0) == [-3.0]
    end
end
