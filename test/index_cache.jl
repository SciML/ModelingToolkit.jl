using ModelingToolkit, SymbolicIndexingInterface, SciMLStructures
using ModelingToolkit: t_nounits as t

# Ensure indexes of array symbolics are cached appropriately
@variables x(t)[1:2]
@named sys = ODESystem(Equation[], t, [x], [])
sys1 = complete(sys)
@named sys = ODESystem(Equation[], t, [x...], [])
sys2 = complete(sys)
for sys in [sys1, sys2]
    for (sym, idx) in [(x, 1:2), (x[1], 1), (x[2], 2)]
        @test is_variable(sys, sym)
        @test variable_index(sys, sym) == idx
    end
end

@variables x(t)[1:2, 1:2]
@named sys = ODESystem(Equation[], t, [x], [])
sys1 = complete(sys)
@named sys = ODESystem(Equation[], t, [x...], [])
sys2 = complete(sys)
for sys in [sys1, sys2]
    @test is_variable(sys, x)
    @test variable_index(sys, x) == [1 3; 2 4]
    for i in eachindex(x)
        @test is_variable(sys, x[i])
        @test variable_index(sys, x[i]) == variable_index(sys, x)[i]
    end
end

# Ensure Symbol to symbolic map is correct
@parameters p1 p2[1:2] p3::String
@variables x(t) y(t)[1:2] z(t)

@named sys = ODESystem(Equation[], t, [x, y, z], [p1, p2, p3])
sys = complete(sys)

ic = ModelingToolkit.get_index_cache(sys)

@test isequal(ic.symbol_to_variable[:p1], p1)
@test isequal(ic.symbol_to_variable[:p2], p2)
@test isequal(ic.symbol_to_variable[:p3], p3)
@test isequal(ic.symbol_to_variable[:x], x)
@test isequal(ic.symbol_to_variable[:y], y)
@test isequal(ic.symbol_to_variable[:z], z)

@testset "tunable_parameters is ordered" begin
    @parameters p q[1:3] r[1:2, 1:2] s [tunable = false]
    @named sys = ODESystem(Equation[], t, [], [p, q, r, s])
    sys = complete(sys)
    @test all(splat(isequal), zip(tunable_parameters(sys), parameters(sys)[1:3]))

    offset = 1
    for par in tunable_parameters(sys)
        idx = parameter_index(sys, par)
        @test idx.portion isa SciMLStructures.Tunable
        if Symbolics.isarraysymbolic(par)
            @test vec(idx.idx) == offset:(offset + length(par) - 1)
        else
            @test idx.idx == offset
        end
        offset += length(par)
    end
end

@testset "reorder_dimension_by_tunables" begin
    @parameters p q[1:3] r[1:2, 1:2] s [tunable = false]
    @named sys = ODESystem(Equation[], t, [], [p, q, r, s])
    src = ones(8)
    dst = zeros(8)
    # system must be complete...
    @test_throws ArgumentError reorder_dimension_by_tunables!(dst, sys, src, [p, q, r])
    @test_throws ArgumentError reorder_dimension_by_tunables(sys, src, [p, q, r])
    sys = complete(sys; split = false)
    # with split = true...
    @test_throws ArgumentError reorder_dimension_by_tunables!(dst, sys, src, [p, q, r])
    @test_throws ArgumentError reorder_dimension_by_tunables(sys, src, [p, q, r])
    sys = complete(sys)
    # and the arrays must have matching size
    @test_throws ArgumentError reorder_dimension_by_tunables!(
        zeros(2, 4), sys, src, [p, q, r])

    ps = MTKParameters(sys, [p => 1.0, q => 3ones(3), r => 4ones(2, 2), s => 0.0])
    src = ps.tunable
    reorder_dimension_by_tunables!(dst, sys, src, [q, r, p])
    @test dst ≈ vcat(3ones(3), 4ones(4), 1.0)
    @test reorder_dimension_by_tunables(sys, src, [r, p, q]) ≈ vcat(4ones(4), 1.0, 3ones(3))
    reorder_dimension_by_tunables!(dst, sys, src, [q[1], r[:, 1], q[2], r[:, 2], q[3], p])
    @test dst ≈ vcat(3.0, 4ones(2), 3.0, 4ones(2), 3.0, 1.0)
    src = stack([copy(ps.tunable) for i in 1:5]; dims = 1)
    dst = zeros(size(src))
    reorder_dimension_by_tunables!(dst, sys, src, [r, q, p]; dim = 2)
    @test dst ≈ stack([vcat(4ones(4), 3ones(3), 1.0) for i in 1:5]; dims = 1)
end

mutable struct ParamTest
    y::Any
end
(pt::ParamTest)(x) = pt.y - x
@testset "Issue#3215: Callable discrete parameter" begin
    function update_affect!(integ, u, p, ctx)
        integ.p[p.p_1].y = integ.t
    end

    tp1 = typeof(ParamTest(1))
    @parameters (p_1::tp1)(..) = ParamTest(1)
    @variables x(ModelingToolkit.t_nounits) = 0

    event1 = [1.0, 2, 3] => (update_affect!, [], [p_1], [p_1], nothing)

    @named sys = ODESystem([
            ModelingToolkit.D_nounits(x) ~ p_1(x)
        ],
        ModelingToolkit.t_nounits;
        discrete_events = [event1]
    )
    ss = @test_nowarn complete(sys)
    @test length(parameters(ss)) == 1
    @test !is_timeseries_parameter(ss, p_1)
end
