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

@testset "Variable discovery from `initial_conditions` and `bindings` kwargs (time-dependent)" begin
    @variables x(t) y(t)   # y is NOT in any equation
    @parameters p q         # q is NOT in any equation

    # Only x and p appear in the equations
    eqs = [D(x) ~ p * x]

    # y appears only as the RHS of an initial_conditions pair
    # q appears only as the RHS of a bindings pair
    # p appears only as the LHS of a bindings pair (already in eqs, but verifying no breakage)
    @named sys = System(eqs, t; initial_conditions = [x => y], bindings = [p => q])

    @test y ∈ Set(unknowns(sys))
    @test q ∈ Set(parameters(sys))
end

@testset "Variable discovery from `initial_conditions` and `bindings` kwargs (time-independent)" begin
    @variables x y          # y is NOT in any equation
    @parameters p q         # q is NOT in any equation

    # Only x and p appear in the equations
    eqs = [0 ~ p * x + 1]

    @named sys = System(eqs; initial_conditions = [x => y], bindings = [p => q])

    @test y ∈ Set(unknowns(sys))
    @test q ∈ Set(parameters(sys))
end
