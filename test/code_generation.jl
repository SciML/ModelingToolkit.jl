using ModelingToolkit, OrdinaryDiffEq, SymbolicIndexingInterface
using ModelingToolkit: t_nounits as t, D_nounits as D

@testset "`generate_custom_function`" begin
    @variables x(t) y(t)[1:3]
    @parameters p1=1.0 p2[1:3]=[1.0, 2.0, 3.0] p3::Int=1 p4::Bool=false

    sys = complete(ODESystem(Equation[], t, [x; y], [p1, p2, p3, p4]; name = :sys))
    u0 = [1.0, 2.0, 3.0, 4.0]
    p = ModelingToolkit.MTKParameters(sys, [])

    fn1 = generate_custom_function(
        sys, x + y[1] + p1 + p2[1] + p3 * t; expression = Val(false))
    @test fn1(u0, p, 0.0) == 5.0

    fn2 = generate_custom_function(
        sys, x + y[1] + p1 + p2[1] + p3 * t, [x], [p1, p2, p3]; expression = Val(false))
    @test fn1(u0, p, 0.0) == 5.0

    fn3_oop, fn3_iip = generate_custom_function(
        sys, [x + y[2], y[3] + p2[2], p1 + p3, 3t]; expression = Val(false))

    buffer = zeros(4)
    fn3_iip(buffer, u0, p, 1.0)
    @test buffer == [4.0, 6.0, 2.0, 3.0]
    @test fn3_oop(u0, p, 1.0) == [4.0, 6.0, 2.0, 3.0]

    fn4 = generate_custom_function(sys, ifelse(p4, p1, p2[2]); expression = Val(false))
    @test fn4(u0, p, 1.0) == 2.0
    fn5 = generate_custom_function(sys, ifelse(!p4, p1, p2[2]); expression = Val(false))
    @test fn5(u0, p, 1.0) == 1.0

    @variables x y[1:3]
    sys = complete(NonlinearSystem(Equation[], [x; y], [p1, p2, p3, p4]; name = :sys))
    p = MTKParameters(sys, [])

    fn1 = generate_custom_function(sys, x + y[1] + p1 + p2[1] + p3; expression = Val(false))
    @test fn1(u0, p) == 6.0

    fn2 = generate_custom_function(
        sys, x + y[1] + p1 + p2[1] + p3, [x], [p1, p2, p3]; expression = Val(false))
    @test fn1(u0, p) == 6.0

    fn3_oop, fn3_iip = generate_custom_function(
        sys, [x + y[2], y[3] + p2[2], p1 + p3]; expression = Val(false))

    buffer = zeros(3)
    fn3_iip(buffer, u0, p)
    @test buffer == [4.0, 6.0, 2.0]
    @test fn3_oop(u0, p, 1.0) == [4.0, 6.0, 2.0]

    fn4 = generate_custom_function(sys, ifelse(p4, p1, p2[2]); expression = Val(false))
    @test fn4(u0, p, 1.0) == 2.0
    fn5 = generate_custom_function(sys, ifelse(!p4, p1, p2[2]); expression = Val(false))
    @test fn5(u0, p, 1.0) == 1.0
end

@testset "Non-standard array variables" begin
    @variables x(t)
    @parameters p[0:2] (f::Function)(..)
    @mtkbuild sys = ODESystem(D(x) ~ p[0] * x + p[1] * t + p[2] + f(p), t)
    prob = ODEProblem(sys, [x => 1.0], (0.0, 1.0), [p => [1.0, 2.0, 3.0], f => sum])
    @test prob.ps[p] == [1.0, 2.0, 3.0]
    @test prob.ps[p[0]] == 1.0
    sol = solve(prob, Tsit5())
    @test SciMLBase.successful_retcode(sol)

    @testset "Array split across buffers" begin
        @variables x(t)[0:2]
        @parameters p[1:2] (f::Function)(..)
        @named sys = ODESystem(
            [D(x[0]) ~ p[1] * x[0] + x[2], D(x[1]) ~ p[2] * f(x) + x[2]], t)
        sys, = structural_simplify(sys, ([x[2]], []))
        @test is_parameter(sys, x[2])
        prob = ODEProblem(sys, [x[0] => 1.0, x[1] => 1.0], (0.0, 1.0),
            [p => ones(2), f => sum, x[2] => 2.0])
        sol = solve(prob, Tsit5())
        @test SciMLBase.successful_retcode(sol)
    end
end

@testset "scalarized array observed calling same function multiple times" begin
    @variables x(t) y(t)[1:2]
    @parameters foo(::Real)[1:2]
    val = Ref(0)
    function _tmp_fn2(x)
        val[] += 1
        return [x, 2x]
    end
    @mtkbuild sys = ODESystem([D(x) ~ y[1] + y[2], y ~ foo(x)], t)
    @test length(equations(sys)) == 1
    @test length(ModelingToolkit.observed(sys)) == 3
    prob = ODEProblem(sys, [x => 1.0], (0.0, 1.0), [foo => _tmp_fn2])
    val[] = 0
    @test_nowarn prob.f(prob.u0, prob.p, 0.0)
    @test val[] == 1

    @testset "CSE in equations(sys)" begin
        val[] = 0
        @variables z(t)[1:2]
        @mtkbuild sys = ODESystem(
            [D(y) ~ foo(x), D(x) ~ sum(y), zeros(2) ~ foo(prod(z))], t)
        @test length(equations(sys)) == 5
        @test length(ModelingToolkit.observed(sys)) == 0
        prob = ODEProblem(
            sys, [y => ones(2), z => 2ones(2), x => 3.0], (0.0, 1.0), [foo => _tmp_fn2])
        val[] = 0
        @test_nowarn prob.f(prob.u0, prob.p, 0.0)
        @test val[] == 2
    end
end
