using ModelingToolkit
using Graphs
using DiffEqBase
using Test
using UnPack
using OrdinaryDiffEq
using LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D

# Define some variables
@parameters L g
@variables x(t) y(t) z(t) w(t) T(t)

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

eqs2 = [D(D(x)) ~ T * x,
    D(D(y)) ~ T * y - g,
    0 ~ x^2 + y^2 - L^2]
pendulum2 = System(eqs2, t, [x, y, T], [L, g], name = :pendulum)

eqs = [D(x) ~ w,
    D(y) ~ z,
    D(w) ~ T * x,
    D(z) ~ T * y - g,
    0 ~ x^2 + y^2 - L^2]
pendulum = System(eqs, t, [x, y, w, z, T], [L, g], name = :pendulum)

let sys = mtkcompile(pendulum2)
    @test length(equations(sys)) == 5
    @test length(unknowns(sys)) == 5

    ivs = [
        x => sqrt(2) / 2,
        y => sqrt(2) / 2,
        L => 1.0,
        g => 9.8
    ]

    prob_auto = ODEProblem(sys, ivs, (0.0, 0.5), guesses = [T => 0.0])
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
    sys = complete(mtkcompile(pend; dummy_derivative = false))
    prob = ODEProblem(
        sys, [x => 1, y => 0, D(x) => 0.0, g => 1], (0.0, 10.0), guesses = [λ => 0.0])
    sol = solve(prob, Rodas5P())
    @test SciMLBase.successful_retcode(sol)
    @test sol[x ^ 2 + y ^ 2][end] < 1.1
end
