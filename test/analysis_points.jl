using ModelingToolkit, ModelingToolkitStandardLibrary
using Test
using ModelingToolkit: t_nounits as t, D_nounits as D

@testset "AnalysisPoint is lowered to `connect`" begin
    @named P = FirstOrder(k = 1, T = 1)
    @named C = Gain(; k = -1)
    
    ap = AnalysisPoint(:plant_input)
    eqs = [connect(P.output, C.input)
           connect(C.output, ap, P.input)]
    sys_ap = ODESystem(eqs, t, systems = [P, C], name = :hej)
    sys_ap2 = @test_nowarn expand_connections(sys)

    @test all(eq -> !(eq.lhs isa AnalysisPoint), equations(sys_ap2))

    eqs = [connect(P.output, C.input)
           connect(C.output, P.input)]
    sys_normal = ODESystem(eqs, t, systems = [P, C], name = :hej)
    sys_normal2 = @test_nowarn expand_connections(sys)

    @test isequal(sys_ap2, sys_normal2)
end
