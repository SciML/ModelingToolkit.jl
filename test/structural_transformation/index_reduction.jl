using ModelingToolkit
using Graphs
using DiffEqBase
using Test
using UnPack
using OrdinaryDiffEq
using LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D

# Simple pendulum in cartesian coordinates
eqs = [D(x) ~ w,
    D(y) ~ z,
    D(w) ~ T * x,
    D(z) ~ T * y - g,
    0 ~ x^2 + y^2 - L^2]
pendulum = System(eqs, t, [x, y, w, z, T], [L, g], name = :pendulum)

state = TearingState(pendulum)
@unpack graph, var_to_diff = state.structure
@test StructuralTransformations.maximal_matching(graph, eq -> true,
    v -> var_to_diff[v] === nothing) ==
      map(x -> x == 0 ? StructuralTransformations.unassigned : x,
    [3, 4, 2, 5, 0, 0, 0, 0, 0])

# Define some variables
@parameters L g
@variables x(t) y(t) T(t)

eqs2 = [D(D(x)) ~ T * x,
    D(D(y)) ~ T * y - g,
    0 ~ x^2 + y^2 - L^2]
pendulum2 = System(eqs2, t, [x, y, T], [L, g], name = :pendulum)

@test_skip begin
    let pss_pendulum2 = partial_state_selection(pendulum2)
        length(equations(pss_pendulum2)) <= 6
    end
end

eqs = [D(x) ~ w,
    D(y) ~ z,
    D(w) ~ T * x,
    D(z) ~ T * y - g,
    0 ~ x^2 + y^2 - L^2]
pendulum = System(eqs, t, [x, y, w, z, T], [L, g], name = :pendulum)

let pss_pendulum = partial_state_selection(pendulum)
    # This currently selects `T` rather than `x` at top level. Needs tearing priorities to fix.
    @test_broken length(equations(pss_pendulum)) == 3
end

let sys = structural_simplify(pendulum2)
    @test length(equations(sys)) == 5
    @test length(unknowns(sys)) == 5

    u0 = [
        x => sqrt(2) / 2,
        y => sqrt(2) / 2
    ]
    p = [
        L => 1.0,
        g => 9.8
    ]

    prob_auto = ODEProblem(sys, u0, (0.0, 0.5), p, guesses = [T => 0.0])
    sol = solve(prob_auto, FBDF())
    @test sol.retcode == ReturnCode.Success
    @test norm(sol[x] .^ 2 + sol[y] .^ 2 .- 1) < 1e-2
end

let
    @parameters g
    @variables x(t) [state_priority = 10] y(t) λ(t)

    eqs = [D(D(x)) ~ λ * x
           D(D(y)) ~ λ * y - g
           x^2 + y^2 ~ 1]
    @named pend = System(eqs, t)
    sys = complete(structural_simplify(pend; dummy_derivative = false))
    prob = ODEProblem(
        sys, [x => 1, y => 0, D(x) => 0.0], (0.0, 10.0), [g => 1], guesses = [λ => 0.0])
    sol = solve(prob, Rodas5P())
    @test SciMLBase.successful_retcode(sol)
    @test sol[x^2 + y^2][end] < 1.1
end
