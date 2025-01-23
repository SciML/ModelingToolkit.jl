using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@testset "Input map validation" begin
    @variables X(t)
    @parameters p d
    eqs = [D(X) ~ p - d*X]
    @mtkbuild osys = ODESystem(eqs, t)
    
    p = "I accidentally renamed p"
    u0 = [X => 1.0]
    ps = [p => 1.0, d => 0.5]
    @test_throws ModelingToolkit.BadKeyError oprob = ODEProblem(osys, u0, (0.0, 1.0), ps)
    
    ps = [p => 1.0, d => 0.5, "Random stuff" => 3.0]
    @test_throws ModelingToolkit.BadKeyError oprob = ODEProblem(osys, u0, (0.0, 1.0), ps)

    u0 = [:X => 1.0, "random" => 3.0]
    @test_throws ModelingToolkit.BadKeyError oprob = ODEProblem(osys, u0, (0.0, 1.0), ps)

    @parameters k
    ps = [p => 1., d => 0.5, k => 3.]
    @test_throws ModelingToolkit.BadKeyError oprob = ODEProblem(osys, u0, (0.0, 1.0), ps)
end
