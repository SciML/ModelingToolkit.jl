using ModelingToolkitBase
using SymbolicIndexingInterface: SymbolicIndexingInterface as SII
using Test
MT = ModelingToolkitBase

@independent_variables t
@variables x
struct MyNLS <: MT.AbstractSystem
    name::Any
    systems::Any
end
tmp = independent_variables(MyNLS("sys", []))
@test tmp == []

struct MyTDS <: MT.AbstractSystem
    iv::Any
    name::Any
    systems::Any
end
iv = independent_variables(MyTDS(t, "sys", []))
@test all(isequal.(iv, [t]))

struct ScalarIVSystem <: MT.AbstractSystem
    time::Any
end
MT.independent_variable(sys::ScalarIVSystem) = getfield(sys, :time)
scalar_iv_sys = ScalarIVSystem(t)
@test isequal(independent_variable(scalar_iv_sys), t)
@test isequal(independent_variables(scalar_iv_sys), [t])

struct MyMVS <: MT.AbstractSystem
    ivs::Any
    name::Any
    systems::Any
end
ivs = independent_variables(MyMVS([t, x], "sys", []))
@test all(isequal.(ivs, [t, x]))

@testset "all_symbols on non-complete system" begin
    @variables y(t)
    @parameters p1 p2 = 2p1
    sys = MT.System(Equation[], t, [y], [p1, p2]; name = :sys)

    # Should not throw on a non-complete system
    syms = SII.all_symbols(sys)
    @test any(isequal(y), syms)
    @test any(isequal(p1), syms)
    @test any(isequal(p2), syms)
    @test any(isequal(t), syms)

    # After completing, bound parameters should also appear
    csys = complete(sys)
    csyms = SII.all_symbols(csys)
    @test any(isequal(csys.y), csyms)
    @test any(isequal(csys.p1), csyms)
    @test any(isequal(t), csyms)

    # p2 is a bound parameter; it should still be in all_symbols
    @test any(isequal(csys.p2), csyms)
    @test any(isequal(csys.p2), collect(MT.bound_parameters(csys)))
end
