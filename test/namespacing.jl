using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D, iscomplete, does_namespacing

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
