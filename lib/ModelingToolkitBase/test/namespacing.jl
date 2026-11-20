using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D, iscomplete, does_namespacing,
    renamespace
using OrdinaryDiffEqTsit5
using Test

@variables x(t)
@parameters p
sys = System(D(x) ~ p * x, t; name = :inner)
@test !iscomplete(sys)
@test does_namespacing(sys)

csys = complete(sys)
@test iscomplete(csys)
@test !does_namespacing(csys)

nsys = toggle_namespacing(sys, false)
@test !iscomplete(nsys)
@test !does_namespacing(nsys)

@test isequal(x, csys.x)
@test isequal(x, nsys.x)
@test !isequal(x, sys.x)
@test isequal(p, csys.p)
@test isequal(p, nsys.p)
@test !isequal(p, sys.p)

@test_throws ["namespacing", "inner"] System(
    Equation[], t; systems = [nsys], name = :a
)

@testset "Vector namespaces" begin
    outer = System(Equation[], t; name = :outer)
    names = [outer, sys]
    for value in (:value, x, D(x) ~ x, [ModelingToolkitBase.unwrap(x)])
        expected = renamespace(outer, renamespace(sys, value))
        @test isequal(renamespace(names, value), expected)
    end
end

@testset "Variables of variables" begin
    @variables x(t) y(x)
    @named inner = System([D(x) ~ x, y ~ 2x + 1], t)
    @test issetequal(unknowns(inner), [x, y])
    ss = mtkcompile(inner)
    @test issetequal(unknowns(ss), [x])

    @named sys = System(Equation[], t; systems = [inner])
    xx, yy = let sys = inner
        xx = renamespace(sys, x)
        yy = only(@variables y(xx))
        xx, renamespace(sys, yy)
    end
    @test issetequal(unknowns(sys), [xx, yy])
    ss = mtkcompile(sys)
    @test isequal(unknowns(ss), [xx])
end

@testset "Array slices" begin
    @variables u(t)[1:5]
    slice = ModelingToolkitBase.unwrap(u)[2:3]
    @named rod = System(Equation[], t)
    @test isequal(renamespace(rod, slice), ModelingToolkitBase.unwrap(renamespace(rod, u))[2:3])

    eqs = [
        D(u[1]) ~ 0,
        D(u[2:4]) ~ u[1:3] .- 2 .* u[2:4] .+ u[3:5],
        D(u[5]) ~ 0,
    ]
    @named inner = System(eqs, t)
    @named outer = System(Equation[], t; systems = [inner])
    u0 = [1.0, 0.0, 0.0, 0.0, 2.0]
    flat = solve(ODEProblem(mtkcompile(inner), [u => u0], (0.0, 1.0)), Tsit5())
    ss = mtkcompile(outer)
    nested = solve(ODEProblem(ss, [ss.inner.u => u0], (0.0, 1.0)), Tsit5())
    @test SciMLBase.successful_retcode(nested)
    @test nested[ss.inner.u][end] ≈ flat[u][end]
end
