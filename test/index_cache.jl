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

