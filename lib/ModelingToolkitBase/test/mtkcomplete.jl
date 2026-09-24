using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D, iscomplete

@testset "@mtkcomplete creates a named, complete system" begin
    @variables x(t)
    @mtkcomplete sys = System([D(x) ~ 2x], t)
    @test nameof(sys) == :sys
    @test iscomplete(sys)
end

@testset "@mtkcomplete with kwargs" begin
    @variables x(t) y(t) z(t)
    @mtkcomplete sys = System([D(x) ~ 2x, z ~ y + x], t; observed = [y ~ 2x + 3])
    @test nameof(sys) == :sys
    @test iscomplete(sys)
    @test isequal(observed(sys), [y ~ 2x + 3])
end

@testset "@mtkcomplete for time-independent system" begin
    @variables x(t)
    @mtkcomplete sys = System([x ~ 1])
    @test nameof(sys) == :sys
    @test iscomplete(sys)
end
