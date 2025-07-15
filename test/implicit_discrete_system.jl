using ModelingToolkit, SymbolicIndexingInterface, Test
using ModelingToolkit: t_nounits as t
using StableRNGs

k = ShiftIndex(t)
rng = StableRNG(22525)

@testset "Correct ImplicitDiscreteFunction" begin
    @variables x(t) = 1
    @mtkcompile sys = System([x(k) ~ x(k) * x(k - 1) - 3], t)
    tspan = (0, 10)

    # u[2] - u_next[1]
    # -3 - u_next[2] + u_next[2]*u_next[1]
    f = ImplicitDiscreteFunction(sys)
    u_next = [3.0, 1.5]
    @test f(u_next, [2.0, 3.0], [], t) ≈ [0.0, 0.0]
    u_next = [0.0, 0.0]
    @test f(u_next, [2.0, 3.0], [], t) ≈ [3.0, -3.0]

    resid = rand(2)
    f(resid, u_next, [2.0, 3.0], [], t)
    @test resid ≈ [3.0, -3.0]

    prob = ImplicitDiscreteProblem(sys, [x(k - 1) => 3.0], tspan)
    @test prob.u0 == [3.0, 1.0]
    prob = ImplicitDiscreteProblem(sys, [], tspan)
    @test prob.u0 == [1.0, 1.0]
    @variables x(t)
    @mtkcompile sys = System([x(k) ~ x(k) * x(k - 1) - 3], t)
    @test_throws ModelingToolkit.MissingGuessError prob=ImplicitDiscreteProblem(
        sys, [], tspan)
end

@testset "System with algebraic equations" begin
    @variables x(t) y(t)
    eqs = [x(k) ~ x(k - 1) + x(k - 2),
        x^2 ~ 1 - y^2]
    @mtkcompile sys = System(eqs, t)
    f = ImplicitDiscreteFunction(sys)

    function correct_f(u_next, u, p, t)
        [u[2] - u_next[1],
            u[1] + u[2] - u_next[2],
            1 - (u_next[1] + u_next[2])^2 - u_next[3]^2]
    end

    reorderer = getu(sys, [x(k - 2), x(k - 1), y])

    for _ in 1:10
        u_next = rand(rng, 3)
        u = rand(rng, 3)
        @test correct_f(u_next, u, [], 0.0) ≈ f(u_next, u, [], 0.0)
    end

    # Initialization is satisfied.
    prob = ImplicitDiscreteProblem(
        sys, [x(k - 1) => 0.3, x(k - 2) => 0.4], (0, 10), guesses = [y => 1])
    @test length(equations(prob.f.initialization_data.initializeprob.f.sys)) == 1
end

@testset "Handle observables in function codegen" begin
    # Observable appears in differential equation
    @variables x(t) y(t) z(t)
    eqs = [x(k) ~ x(k - 1) + x(k - 2),
        y(k) ~ x(k) + x(k - 2) * z(k - 1),
        x + y + z ~ 2]
    @mtkcompile sys = System(eqs, t)
    @test length(unknowns(sys)) == length(equations(sys)) == 3
    @test occursin(
        "var\"y(t)\"", string(ImplicitDiscreteFunction(sys; expression = Val{true})))

    # Shifted observable that appears in algebraic equation is properly handled.
    eqs = [z(k) ~ x(k) + sin(x(k)),
        y(k) ~ x(k - 1) + x(k - 2),
        z(k) * x(k) ~ 3]
    @mtkcompile sys = System(eqs, t)
    @test occursin("var\"Shift(t, 1)(z(t))\"",
        string(ImplicitDiscreteFunction(sys; expression = Val{true})))
end
