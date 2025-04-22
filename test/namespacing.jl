using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D, iscomplete, does_namespacing

@testset "System" begin
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
        Equation[], t; systems = [nsys], name = :a)
end

@testset "SDESystem" begin
    @variables x(t)
    @parameters p
    sys = SDESystem([D(x) ~ p * x], [x], t, [x], [p]; name = :inner)
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

    @test_throws ["namespacing", "inner"] SDESystem(
        Equation[], [], t; systems = [nsys], name = :a)
end

@testset "DiscreteSystem" begin
    @variables x(t)
    @parameters p
    k = ShiftIndex(t)
    sys = DiscreteSystem([x(k) ~ p * x(k - 1)], t; name = :inner)
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

    @test_throws ["namespacing", "inner"] DiscreteSystem(
        Equation[], t; systems = [nsys], name = :a)
end

@testset "ImplicitDiscreteSystem" begin
    @variables x(t)
    @parameters p
    k = ShiftIndex(t)
    sys = System([x(k) ~ p + x(k - 1) * x(k)], t; name = :inner)
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
        Equation[], t; systems = [nsys], name = :a)
end

@testset "NonlinearSystem" begin
    @variables x
    @parameters p
    sys = System([x ~ p * x^2 + 1]; name = :inner)
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
        Equation[]; systems = [nsys], name = :a)
end

@testset "OptimizationSystem" begin
    @variables x
    @parameters p
    sys = OptimizationSystem(p * x^2 + 1; name = :inner)
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

    @test_throws ["namespacing", "inner"] OptimizationSystem(
        []; systems = [nsys], name = :a)
end

@testset "ConstraintsSystem" begin
    @variables x
    @parameters p
    sys = ConstraintsSystem([x^2 + p ~ 0], [x], [p]; name = :inner)
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

    @test_throws ["namespacing", "inner"] ConstraintsSystem(
        [], [], []; systems = [nsys], name = :a)
end
