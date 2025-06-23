using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@testset "Input map validation" begin
    import ModelingToolkit: InvalidKeyError, MissingParametersError
    @variables X(t)
    @parameters p d
    eqs = [D(X) ~ p - d * X]
    @mtkcompile osys = System(eqs, t)

    p = "I accidentally renamed p"
    u0 = [X => 1.0]
    ps = [p => 1.0, d => 0.5]
    @test_throws MissingParametersError oprob=ODEProblem(osys, u0, (0.0, 1.0), ps)

    @parameters p d
    ps = [p => 1.0, d => 0.5, "Random stuff" => 3.0]
    @test_throws InvalidKeyError oprob=ODEProblem(osys, u0, (0.0, 1.0), ps)

    u0 = [:X => 1.0, "random" => 3.0]
    @test_throws InvalidKeyError oprob=ODEProblem(osys, u0, (0.0, 1.0), ps)

    @variables x(t) y(t) z(t)
    @parameters a b c d
    eqs = [D(x) ~ x * a, D(y) ~ y * c, D(z) ~ b + d]
    @mtkcompile sys = System(eqs, t)
    pmap = [a => 1, b => 2, c => 3, d => 4, "b" => 2]
    u0map = [x => 1, y => 2, z => 3]
    @test_throws InvalidKeyError ODEProblem(sys, u0map, (0.0, 1.0), pmap)

    pmap = [a => 1, b => 2, c => 3, d => 4]
    u0map = [x => 1, y => 2, z => 3, :0 => 3]
    @test_throws InvalidKeyError ODEProblem(sys, u0map, (0.0, 1.0), pmap)
end
