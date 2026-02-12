using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D

@testset "User-provided `observed` is respected" begin
    @variables x(t) y(t) z(t)
    @mtkcompile sys = System([D(x) ~ 2x, z ~ y + x], t; observed = [y ~ 2x + 3])
    @test length(equations(sys)) == 1
    @test isequal(observed(sys), [y ~ 2x + 3, z ~ y + x])
end
