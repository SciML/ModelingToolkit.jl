using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D, SymbolicContinuousCallback, SymbolicDiscreteCallback
using Symbolics: SConst

@testset "Simple metadata bindings" begin
    @variables x(t) = 1 y(t) = x
    @parameters p = 1 q = p
    @named sys = System([D(x) ~ p * x, D(y) ~ q * y], t)
    @test isequal(bindings(sys), Dict(y => x, q => p))
    @test isequal(initial_conditions(sys), Dict(x => SConst(1), p => SConst(1)))
end

@testset "Array bindings" begin
    @variables x(t)[1:2] = [1.0, 2.0] y(t)[1:2] = [x[1], 2.0]
    @parameters p[1:2] = [2.0, 3.0] q[1:2] = [p[1], 4.0]
    @named sys = System(Equation[], t, [x, y], [p, q])
    @test isequal(bindings(sys), Dict(y => SConst([x[1], 2.0]), q => SConst([p[1], 4.0])))
    @test isequal(initial_conditions(sys), Dict(x => SConst([1.0, 2.0]), p => SConst([2.0, 3.0])))
end

@testset "Prefer keyword over metadata" begin
    @variables x(t) = 1 y(t) = x
    @parameters p = 1 q = p
    @named sys = System([D(x) ~ p * x, D(y) ~ q * y], t;
                        bindings = [x => y, p => q, y => nothing, q => nothing],
                        initial_conditions = [y => 2, q => 2, x => nothing, p => nothing])
    @test isequal(bindings(sys), Dict(x => y, p => q))
    @test isequal(initial_conditions(sys), Dict(y => SConst(2), q => SConst(2)))
end

@testset "Arrays are atomic" begin
    @variables x(t)[1:2]
    @test_throws ModelingToolkit.IndexedArrayKeyError System(D(x) ~ x, t; bindings = [x[1] => x[2]], name = :a)
    @test_throws ModelingToolkit.IndexedArrayKeyError System(D(x) ~ x, t; initial_conditions = [x[1] => 1], name = :a)
    @test_throws ModelingToolkit.IndexedArrayKeyError System(D(x) ~ x, t; guesses = [x[1] => 1], name = :a)
end

@testset "Parameter bindings cannot involve variables" begin
    @variables x(t)
    @parameters p = 2x + 1
    @test_throws ["parameter p", "encountered binding", "non-parameter"] System(D(x) ~ x, t, [x], [p]; name = :a)
end

@testset "Discrete unknowns eventually become parameters" begin
    @variables x(t)
    @discretes d1(t) = x d2(t) = d1
    cev = SymbolicContinuousCallback([x ~ 1.0], [d1 ~ Pre(x)]; discrete_parameters = [d1])
    dev = SymbolicDiscreteCallback(1.0, [d2 ~ Pre(x)]; discrete_parameters = [d2])
    @named sys = System([D(x) ~ x], t, [x, d1, d2], []; continuous_events = [cev], discrete_events = [dev])
    @test d1 in Set(unknowns(sys))
    @test d2 in Set(unknowns(sys))
    csys = ModelingToolkit.discrete_unknowns_to_parameters(sys)
    @test d1 in Set(parameters(csys))
    @test d2 in Set(parameters(csys))
    ts = TearingState(sys)
    ssys = ts.sys
    @test d1 in Set(parameters(ssys))
    @test d2 in Set(parameters(ssys))
end

@testset "`missing` bindings are bindings and not initial conditions" begin
    @variables x(t)
    @parameters p = missing
    @named sys = System(D(x) ~ x * p, t)
    @test bindings(sys)[p] === ModelingToolkit.COMMON_MISSING
end
