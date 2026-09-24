using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D

@testset "User-provided `observed` is respected" begin
    @variables x(t) y(t) z(t)
    @mtkcompile sys = System([D(x) ~ 2x, z ~ y + x], t; observed = [y ~ 2x + 3])
    @test length(equations(sys)) == 1
    @test isequal(observed(sys), [y ~ 2x + 3, z ~ y + x])

    @testset "In nonlinear system" begin
        @variables X1(t) X2(t) X3(t)
        @parameters Γ[1:1] k1 k2 k3 k4
        nleqs = [
            0 ~ -k1 * X1 + k2 * X2 - k3 * X2 * X1 + (1 // 2) * k4 * (X3^2),
            0 ~ k1 * X1 - k2 * X2 - k3 * X2 * X1 + (1 // 2) * k4 * (X3^2),
        ]
        @mtkcompile sys = System(nleqs; observed = [X3 ~ Γ[1] - X1 - X2])
        @test issetequal(equations(sys), nleqs)
        @test isequal(only(observed(sys)), X3 ~ Γ[1] - X1 - X2)
    end
end
