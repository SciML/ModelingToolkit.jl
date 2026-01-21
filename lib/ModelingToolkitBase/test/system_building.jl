using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D

@testset "`state_priorities` and `irreducibles` kwargs take priority over metadata" begin
    @variables x(t) [state_priority = 3] y(t) [irreducible = false]
    @named sys = System([D(x) ~ 2x + y], t; state_priorities = [x => -4], irreducibles = [y])
    @test state_priorities(sys)[x] == -4
    @test y in irreducibles(sys)
end

@testset "Issue#4138: expressions of `Initial` allowed as parameter defaults" begin
    @variables X1(t) X2(t)
    @parameters k1 k2 Γ = Initial(X1) + Initial(X2)
    eqs = [
        D(X1) ~ -k1 * X1 + k2 * (-X1 + Γ),
        X2 ~ -X1 + Γ,
    ]
    @mtkcompile sys = System(eqs, t)
    @test isequal(bindings(sys)[Γ], Initial(X1) + Initial(X2))
end
