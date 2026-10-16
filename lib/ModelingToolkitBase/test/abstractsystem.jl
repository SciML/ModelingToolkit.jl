using ModelingToolkitBase
using SymbolicIndexingInterface: SymbolicIndexingInterface as SII
import REPL
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

using ModelingToolkitBase: t_nounits as t, D_nounits as D

struct NotASystem <: ModelingToolkitBase.AbstractSystem end

@testset "`extend` rejects systems of different types" begin
    @variables x(t)
    @named sys = System([D(x) ~ x], t)
    @test_throws ArgumentError extend(
        NotASystem(), sys; name = :ext, description = "", gui_metadata = nothing
    )
end

@testset "`extend` keeps constraints" begin
    @variables x(t) y(t)
    @parameters p
    @named sys1 = System([D(x) ~ -x], t; constraints = [x ~ p])
    @named sys2 = System([y ~ 2x], t)
    @test issetequal(MT.get_constraints(extend(sys2, sys1)), [x ~ p])
    @test issetequal(MT.get_constraints(extend(sys1, sys2)), [x ~ p])
end

@testset "`propertynames` skips unknowns and observed without names" begin
    # https://github.com/SciML/ModelingToolkit.jl/issues/923
    @variables x(t)
    sys = System([D(x) ~ x + 1], t, [x, D(x)], []; name = :sys)
    for s in (sys, complete(sys))
        names = propertynames(s)
        @test :x in names && :t in names
        for n in names
            @test_nowarn getproperty(s, n)
        end
    end

    sys_obs = System([D(x) ~ x], t; observed = [D(x) ~ 2x], name = :sys_obs)
    names = propertynames(sys_obs)
    @test :x in names && :t in names
    for n in names
        @test_nowarn getproperty(sys_obs, n)
    end

    # REPL tab completion calls `propertynames` on the value before the dot.
    Core.eval(@__MODULE__, :(issue923_sys = $sys))
    completions = REPL.REPLCompletions.completions(
        "issue923_sys.", length("issue923_sys."), @__MODULE__
    )[1]
    @test "x" in REPL.REPLCompletions.completion_text.(completions)
end
