using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D, iscomplete, does_namespacing,
    renamespace
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
