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
eqs = [
    D(x) ~ w,
    D(y) ~ z,
    D(w) ~ T * x,
    D(z) ~ T * y - g,
    0 ~ x^2 + y^2 - L^2,
]
pendulum = System(eqs, t, [x, y, w, z, T], [L, g], name = :pendulum)

state = TearingState(pendulum)
@unpack graph, var_to_diff = state.structure
@test StructuralTransformations.maximal_matching(
    graph, eq -> true,
    v -> var_to_diff[v] === nothing
) == map(state.fullvars) do v
    if operation(v) isa Differential
        return findfirst(eq -> isequal(eq.lhs, v), equations(state))
    end
    return StructuralTransformations.unassigned
end
eqs2 = [
    D(D(x)) ~ T * x,
    D(D(y)) ~ T * y - g,
    0 ~ x^2 + y^2 - L^2,
]
pendulum2 = System(eqs2, t, [x, y, T], [L, g], name = :pendulum)

eqs = [
    D(x) ~ w,
    D(y) ~ z,
    D(w) ~ T * x,
    D(z) ~ T * y - g,
    0 ~ x^2 + y^2 - L^2,
]
pendulum = System(eqs, t, [x, y, w, z, T], [L, g], name = :pendulum)

let sys = mtkcompile(pendulum2)
    ivs = [
        x => sqrt(2) / 2,
        y => sqrt(2) / 2,
        L => 1.0,
        g => 9.8,
    ]

    prob_auto = ODEProblem(sys, ivs, (0.0, 0.5), guesses = [T => 0.0])
    sol = solve(prob_auto, FBDF())
    @test sol.retcode == ReturnCode.Success
    @test norm(sol[x] .^ 2 + sol[y] .^ 2 .- 1) < 1.0e-2
end

let
    @parameters g
    @variables x(t) [state_priority = 10] y(t) λ(t)

    eqs = [
        D(D(x)) ~ λ * x
        D(D(y)) ~ λ * y - g
        x^2 + y^2 ~ 1
    ]
    @named pend = System(eqs, t)
    sys = complete(mtkcompile(pend; dummy_derivative = false))
    prob = ODEProblem(
        sys, [x => 1, y => 0, D(x) => 0.0, g => 1], (0.0, 10.0), guesses = [λ => 0.0]
    )
    sol = solve(prob, Rodas5P())
    @test SciMLBase.successful_retcode(sol)
    @test sol[x^2 + y^2][end] < 1.1
end

# Index reduction through a struct-typed unknown. `pendulum2` above is the same system
# written with three free-standing scalars, so it is the control: the record version must
# reduce to the same structure and integrate to the same trajectory.
let
    struct PendulumState
        x::Float64
        y::Float64
        λ::Float64
    end
    @symstruct PendulumState

    @parameters L g
    @variables p(t)::PendulumState

    eqs = [D(D(p.x)) ~ p.λ * p.x
           D(D(p.y)) ~ p.λ * p.y - g
           0 ~ p.x^2 + p.y^2 - L^2]
    struct_sys = mtkcompile(System(eqs, t, [p], [L, g]; name = :pendulum_struct))
    scalar_sys = mtkcompile(pendulum2)

    # Same reduction: index 3 -> 1, same number of states and equations.
    @test length(unknowns(struct_sys)) == length(unknowns(scalar_sys))
    @test length(equations(struct_sys)) == length(equations(scalar_sys))
    # Every state is a field access, i.e. a projection of a record rather than a
    # free-standing variable.
    @test all(ModelingToolkitBase.is_symstruct_field, unknowns(struct_sys))

    ivs = [p.x => sqrt(2) / 2, p.y => sqrt(2) / 2, L => 1.0, g => 9.8]
    sivs = [x => sqrt(2) / 2, y => sqrt(2) / 2, L => 1.0, g => 9.8]
    struct_prob = ODEProblem(struct_sys, ivs, (0.0, 0.5); guesses = [p.λ => 0.0])
    scalar_prob = ODEProblem(scalar_sys, sivs, (0.0, 0.5); guesses = [T => 0.0])

    # A surviving constraint means a singular mass matrix in both cases.
    @test struct_prob.f.mass_matrix == scalar_prob.f.mass_matrix
    @test !(struct_prob.f.mass_matrix === I)

    kw = (; saveat = 0.05, reltol = 1e-10, abstol = 1e-10)
    ssol = solve(struct_prob, FBDF(); kw...)
    csol = solve(scalar_prob, FBDF(); kw...)
    @test SciMLBase.successful_retcode(ssol)
    @test ssol.t == csol.t
    @test maximum(abs.(ssol[p.x] .- csol[x])) < 1e-8
    @test maximum(abs.(ssol[p.y] .- csol[y])) < 1e-8
    # The constraint is respected along the whole trajectory.
    @test norm(ssol[p.x] .^ 2 + ssol[p.y] .^ 2 .- 1) < 1.0e-2
end
