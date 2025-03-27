using ModelingToolkit
using Graphs
using DiffEqBase
using Test
using UnPack
using ModelingToolkit: t_nounits as t, D_nounits as D

# Define some variables
@parameters L g
@variables x(t) y(t) w(t) z(t) T(t) xˍt(t) yˍt(t) xˍˍt(t) yˍˍt(t)

eqs2 = [D(D(x)) ~ T * x,
    D(D(y)) ~ T * y - g,
    0 ~ x^2 + y^2 - L^2]
pendulum2 = ODESystem(eqs2, t, [x, y, T], [L, g], name = :pendulum)
lowered_sys = ModelingToolkit.ode_order_lowering(pendulum2)

lowered_eqs = [D(xˍt) ~ T * x,
    D(yˍt) ~ T * y - g,
    D(x) ~ xˍt,
    D(y) ~ yˍt,
    0 ~ x^2 + y^2 - L^2]
@test ODESystem(lowered_eqs, t, [xˍt, yˍt, x, y, T], [L, g], name = :pendulum) ==
      lowered_sys
@test isequal(equations(lowered_sys), lowered_eqs)

# Simple pendulum in cartesian coordinates
eqs = [D(x) ~ w,
    D(y) ~ z,
    D(w) ~ T * x,
    D(z) ~ T * y - g,
    0 ~ x^2 + y^2 - L^2]
pendulum = ODESystem(eqs, t, [x, y, w, z, T], [L, g], name = :pendulum)

state = TearingState(pendulum)
@unpack graph, var_to_diff = state.structure
@test StructuralTransformations.maximal_matching(graph, eq -> true,
    v -> var_to_diff[v] === nothing) ==
      map(x -> x == 0 ? StructuralTransformations.unassigned : x,
    [1, 2, 3, 4, 0, 0, 0, 0, 0])

using ModelingToolkit
@parameters L g
@variables x(t) y(t) w(t) z(t) T(t) xˍt(t) yˍt(t)
idx1_pendulum = [D(x) ~ w,
    D(y) ~ z,
    #0 ~ x^2 + y^2 - L^2,
    D(w) ~ T * x,
    D(z) ~ T * y - g,
    # intermediate 1: 0 ~ 2x*D(x) + 2y*D(y) - 0,
    # intermediate 2(a): 0 ~ 2x*w + 2y*z - 0, (substitute D(x) and D(y))
    #0 ~ 2x*w + 2y*z,
    # D(D(x)) ~ D(w) and substitute the rhs
    D(xˍt) ~ T * x,
    # D(D(y)) ~ D(z) and substitute the rhs
    D(yˍt) ~ T * y - g,
    # 2x*D(D(x)) + 2*D(x)*D(x) + 2y*D(D(y)) + 2*D(y)*D(y) and
    # substitute the rhs
    0 ~ 2x * (T * x) + 2 * xˍt * xˍt + 2y * (T * y - g) + 2 * yˍt * yˍt]
@named idx1_pendulum = ODESystem(idx1_pendulum, t, [x, y, w, z, xˍt, yˍt, T], [L, g])
first_order_idx1_pendulum = complete(ode_order_lowering(idx1_pendulum))

using OrdinaryDiffEq
using LinearAlgebra
prob = ODEProblem(first_order_idx1_pendulum,
    #  [x, y, w, z, xˍt, yˍt, T]
    [1, 0, 0, 0, 0, 0, 0.0],# 0, 0, 0, 0],
    (0, 10.0),
    [1, 9.8])
sol = solve(prob, Rodas5());
#plot(sol, idxs=(1, 2))

new_sys = complete(dae_index_lowering(ModelingToolkit.ode_order_lowering(pendulum2)))

prob_auto = ODEProblem(new_sys,
    [D(x) => 0,
        D(y) => 0,
        x => 1,
        y => 0,
        T => 0.0],
    (0, 100.0),
    [L => 1, g => 9.8])
sol = solve(prob_auto, Rodas5());
#plot(sol, idxs=(x, y))

# Define some variables
@parameters L g
@variables x(t) y(t) T(t)

eqs2 = [D(D(x)) ~ T * x,
    D(D(y)) ~ T * y - g,
    0 ~ x^2 + y^2 - L^2]
pendulum2 = ODESystem(eqs2, t, [x, y, T], [L, g], name = :pendulum)

# Turn into a first order differential equation system
first_order_sys = ModelingToolkit.ode_order_lowering(pendulum2)

# Perform index reduction to get an Index 1 DAE
new_sys = complete(dae_index_lowering(first_order_sys))

u0 = [
    D(x) => 0.0,
    D(y) => 0.0,
    x => 1.0,
    y => 0.0,
    T => 0.0
]

p = [
    L => 1.0,
    g => 9.8
]

prob_auto = ODEProblem(new_sys, u0, (0.0, 10.0), p)
sol = solve(prob_auto, Rodas5());
#plot(sol, idxs=(D(x), y))

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
pendulum = ODESystem(eqs, t, [x, y, w, z, T], [L, g], name = :pendulum)

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
    @named pend = ODESystem(eqs, t)
    sys = complete(structural_simplify(pend; dummy_derivative = false))
    prob = ODEProblem(
        sys, [x => 1, y => 0, D(x) => 0.0], (0.0, 10.0), [g => 1], guesses = [λ => 0.0])
    sol = solve(prob, Rodas5P())
    @test SciMLBase.successful_retcode(sol)
    @test sol[x^2 + y^2][end] < 1.1
end
