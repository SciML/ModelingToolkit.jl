using ModelingToolkit, OrdinaryDiffEq, Test, NonlinearSolve, LinearAlgebra
using BipartiteGraphs
using Symbolics
using OrdinaryDiffEqRosenbrock
using ModelingToolkit: topsort_equations, t_nounits as t, D_nounits as D, unwrap

@variables x(t) y(t) z(t) k(t)
eqs = [
    x ~ y + z
    z ~ 2
    y ~ 2z + k
]
@named sys = System(Equation[], t, [x, y, z, k], [])

sorted_eq = topsort_equations(sys, eqs, unwrap.([x, y, z, k]))

ref_eq = [
    z ~ 2
    y ~ 2z + k
    x ~ y + z
]
@test ref_eq == sorted_eq

@test_throws ArgumentError topsort_equations(
    sys,
    [
        x ~ y + z
        z ~ 2
        y ~ 2z + x
    ], unwrap.([x, y, z, k])
)

@parameters σ ρ β
@variables x(t) y(t) z(t) a(t) u(t) F(t)

test_equal(a, b) = @test isequal(a, b) || isequal(simplify(a), simplify(b))

eqs = [
    D(x) ~ σ * (y - x)
    D(y) ~ x * (ρ - z) - y + β
    0 ~ z - x + y
    0 ~ a + z
    u ~ z + a
]

lorenz1 = System(eqs, t, name = :lorenz1)

lorenz1_aliased = mtkcompile(lorenz1)
io = IOBuffer();
show(io, MIME("text/plain"), lorenz1_aliased);
str = String(take!(io));
@test all(s -> occursin(s, str), ["lorenz1", "Unknowns (2)", "Parameters (3)"])
reduced_eqs = [
    D(x) ~ σ * (y - x)
    D(y) ~ β + (ρ - z) * x - y
]
#test_equal.(equations(lorenz1_aliased), reduced_eqs)
@test isempty(setdiff(unknowns(lorenz1_aliased), [x, y, z]))
#test_equal.(observed(lorenz1_aliased), [u ~ 0
#                                        z ~ x - y
#                                        a ~ -z])

# Multi-System Reduction

@variables s(t)
eqs1 = [
    D(x) ~ σ * (y - x) + F,
    D(y) ~ x * (ρ - z) - u,
    D(z) ~ x * y - β * z,
    u ~ x + y - z,
]

lorenz = name -> System(eqs1, t, name = name)
lorenz1 = lorenz(:lorenz1)
state = TearingState(lorenz1)
@test isempty(setdiff(state.fullvars, [D(x), F, y, x, D(y), u, z, D(z)]))
lorenz2 = lorenz(:lorenz2)

@named connected = System(
    [
        s ~ a + lorenz1.x
        lorenz2.y ~ s
        lorenz1.u ~ lorenz2.F
        lorenz2.u ~ lorenz1.F
    ],
    t, systems = [lorenz1, lorenz2]
)
@test length(Base.propertynames(connected)) == 10 + 1 # + 1 for independent variable
@test isequal((@nonamespace connected.lorenz1.x), x)
__x = x
@unpack lorenz1 = connected
@unpack x = lorenz1
@test isequal(x, __x)

# Reduced Flattened System

reduced_system = mtkcompile(connected)
reduced_system2 = mtkcompile(tearing_substitution(mtkcompile(tearing_substitution(mtkcompile(connected)))))

@test isempty(setdiff(unknowns(reduced_system), unknowns(reduced_system2)))
@test isequal(equations(tearing_substitution(reduced_system)), equations(reduced_system2))
@test isequal(observed(reduced_system), observed(reduced_system2))
@test setdiff(
    unknowns(reduced_system),
    [
        s
        a
        lorenz1.x
        lorenz1.y
        lorenz1.z
        lorenz1.u
        lorenz2.x
        lorenz2.y
        lorenz2.z
        lorenz2.u
    ]
) |> isempty

@test setdiff(
    parameters(reduced_system),
    [
        lorenz1.σ
        lorenz1.ρ
        lorenz1.β
        lorenz2.σ
        lorenz2.ρ
        lorenz2.β
    ]
) |> isempty

@test length(equations(reduced_system)) == 6

pp = [
    lorenz1.σ => 10
    lorenz1.ρ => 28
    lorenz1.β => 8 / 3
    lorenz2.σ => 10
    lorenz2.ρ => 28
    lorenz2.β => 8 / 3
]
u0 = [
    lorenz1.x => 1.0
    lorenz1.y => 0.0
    lorenz1.z => 0.0
    s => 0.0
    lorenz2.x => 1.0
    lorenz2.y => 0.0
    lorenz2.z => 0.0
]
prob1 = ODEProblem(reduced_system, [u0; pp], (0.0, 100.0))
solve(prob1, Rodas5())

prob2 = SteadyStateProblem(reduced_system, [u0; pp])
@test prob2.f.observed(lorenz2.u, prob2.u0, prob2.p) === 1.0

# issue #724 and #716
let
    @variables x(t) u(t) y(t)
    @parameters a b c d
    ol = System([D(x) ~ a * x + b * u; y ~ c * x + d * u], t, name = :ol)
    @variables u_c(t) y_c(t)
    @parameters k_P
    pc = System(Equation[u_c ~ k_P * y_c], t, name = :pc)
    connections = [
        pc.u_c ~ ol.u
        pc.y_c ~ ol.y
    ]
    @named connected = System(connections, t, systems = [ol, pc])
    @test equations(connected) isa Vector{Equation}
    reduced_sys = mtkcompile(connected)
    ref_eqs = [
        D(ol.x) ~ ol.a * ol.x + ol.b * ol.u
        0 ~ pc.k_P * ol.y - ol.u
    ]
    #@test ref_eqs == equations(reduced_sys)
end

# issue #889
let
    @variables x(t)
    @named sys = System([0 ~ D(x) + x], t, [x], [])
    #@test_throws ModelingToolkit.InvalidSystemException ODEProblem(sys, [1.0], (0, 10.0))
    sys = mtkcompile(sys)
    #@test_nowarn ODEProblem(sys, [1.0], (0, 10.0))
end

# NonlinearSystem
@variables u1(t) u2(t) u3(t)
@parameters p
eqs = [
    u1 ~ u2
    u3 ~ u1 + u2 + p
    u3 ~ hypot(u1, u2) * p
]
@named sys = System(eqs, [u1, u2, u3], [p])
reducedsys = mtkcompile(sys)
@test length(observed(reducedsys)) == 2

u0 = [u2 => 1]
pp = [2]
nlprob = NonlinearProblem(reducedsys, [u0; [p => pp[1]]])
reducedsol = solve(nlprob, NewtonRaphson())
residual = fill(100.0, length(unknowns(reducedsys)))
nlprob.f(residual, reducedsol.u, pp)
@test hypot(
    nlprob.f.observed(u2, reducedsol.u, pp),
    nlprob.f.observed(u1, reducedsol.u, pp)
) *
    pp[1] ≈ nlprob.f.observed(u3, reducedsol.u, pp) atol = 1.0e-9

@test all(x -> abs(x) < 1.0e-5, residual)

N = 5
@variables xs[1:N]
A = reshape(1:(N^2), N, N)
eqs = xs ~ A * xs
@named sys′ = System(eqs, [xs], [])
sys = mtkcompile(sys′)
@test isempty(equations(sys)) && length(observed(sys)) == 6 # 5 + 1 for change_origin

# issue 958
@parameters k₁ k₂ k₋₁ E₀
@variables E(t) C(t) S(t) P(t)

eqs = [
    D(E) ~ k₋₁ * C - k₁ * E * S
    D(C) ~ k₁ * E * S - k₋₁ * C - k₂ * C
    D(S) ~ k₋₁ * C - k₁ * E * S
    D(P) ~ k₂ * C
    E₀ ~ E + C
]

@named sys = System(eqs, t, [E, C, S, P], [k₁, k₂, k₋₁, E₀])
@test_throws ModelingToolkit.StateSelection.ExtraEquationsSystemException mtkcompile(sys)

# Example 5 from Pantelides' original paper
params = collect(@parameters y1 y2)
sts = collect(@variables x(t) u1(t) u2(t))
eqs = [
    0 ~ x + sin(u1 + u2)
    D(x) ~ x + y1
    cos(x) ~ sin(y2)
]
@named sys = System(eqs, t, sts, params)
@test_throws ModelingToolkit.StateSelection.InvalidSystemException mtkcompile(sys)

# issue #963
@variables v47(t) v57(t) v66(t) v25(t) i74(t) i75(t) i64(t) i71(t) v1(t) v2(t)

eq = [
    v47 ~ v1
    v47 ~ sin(10t)
    v57 ~ v1 - v2
    v57 ~ 10.0i64
    v66 ~ v2
    v66 ~ 5.0i74
    v25 ~ v2
    i75 ~ 0.005 * D(v25)
    0 ~ i74 + i75 - i64
    0 ~ i64 + i71
]

@named sys0 = System(eq, t)
sys = mtkcompile(sys0)
@test length(equations(sys)) == 1
eq = equations(tearing_substitution(sys))[1]
vv = only(unknowns(sys))
@test isequal(eq.lhs, D(vv))
dvv = ModelingToolkit.value(ModelingToolkit.derivative(eq.rhs, vv))
@test dvv ≈ -60

# Don't reduce inputs
@parameters σ ρ β
@variables x(t) y(t) z(t) [input = true] a(t) u(t) F(t)

eqs = [
    D(x) ~ σ * (y - x)
    D(y) ~ x * (ρ - z) - y + β
    0 ~ a + z
    u ~ z + a
]

lorenz1 = System(eqs, t, name = :lorenz1)
lorenz1_reduced = mtkcompile(lorenz1, inputs = [z], outputs = [])
@test z in Set(parameters(lorenz1_reduced))

# #2064
vars = @variables x(t) y(t) z(t)
eqs = [
    D(x) ~ x
    D(y) ~ y
    D(z) ~ t
]
@named model = System(eqs, t)
sys = mtkcompile(model)
Js = ModelingToolkit.jacobian_sparsity(sys)
@test size(Js) == (3, 3)
@test Js == Diagonal([0, 1, 1])

# MWE for #1722
vars = @variables a(t) w(t) phi(t)
eqs = [
    a ~ D(w)
    w ~ D(phi)
    w ~ sin(t)
]
@named sys = System(eqs, t, vars, [])
ss = alias_elimination(sys)
@test isempty(observed(ss))

@variables x(t) y(t)
@named sys = System(
    [
        D(x) ~ 1 - x,
        D(y) + D(x) ~ 0,
    ], t
)
new_sys = alias_elimination(sys)
@test isempty(observed(new_sys))

@named sys = System(
    [
        D(x) ~ x,
        D(y) + D(x) ~ 0,
    ], t
)
new_sys = alias_elimination(sys)
@test isempty(observed(new_sys))

@named sys = System(
    [
        D(x) ~ 1 - x,
        y + D(x) ~ 0,
    ], t
)
new_sys = alias_elimination(sys)
@test isempty(observed(new_sys))

eqs = [
    x ~ 0
    D(x) ~ x + y
]
@named sys = System(eqs, t, [x, y], [])
ss = mtkcompile(sys)
@test isempty(equations(ss))
dx = ModelingToolkit.default_toterm(unwrap(D(x)))
@test issetequal(observed(ss), [x ~ 0, dx ~ 0, y ~ dx - x])

eqs = [D(D(x)) ~ -x]
@named sys = System(eqs, t, [x], [])
ss = alias_elimination(sys)
@test length(equations(ss)) == length(unknowns(ss)) == 1
ss = mtkcompile(sys)
@test length(equations(ss)) == length(unknowns(ss)) == 2

@testset "Aliases of differential variables with higher state priority are swapped" begin
    @variables x(t) y(t)
    @mtkcompile sys = System([D(x) ~ 2x, y ~ x], t; state_priorities = [y => 10])
    @test isequal(only(unknowns(sys)), y)
end

@testset "positive-priority variable is eliminated when target has higher priority" begin
    # Both x and y have positive state_priority, but y's is higher.
    # Old behaviour: both sticky → neither eliminated.
    # New behaviour: x (lower priority) is eliminated as observed in favour of y.
    @variables x(t) y(t)
    @named sys = System([D(y) ~ -y, x ~ y], t; state_priorities = [x => 1, y => 2])
    state = TearingState(sys)
    ModelingToolkit.eliminate_perfect_aliases!(state)
    @test any(eq -> isequal(eq.lhs, unwrap(x)), state.additional_observed)
    @test any(isequal(unwrap(y)), state.fullvars)
    @test !any(isequal(unwrap(x)), state.fullvars)
end

@testset "warning is emitted when state_priority is tied across alias group" begin
    # x appears in two equations (D(x) ~ -x and x ~ y), y in one (x ~ y),
    # so the equation-count heuristic picks x as target — but the priority tie
    # (both at 5) must still emit a warning regardless of how the tie is broken.
    @variables x(t) y(t)
    @named sys = System([D(x) ~ -x, x ~ y], t; state_priorities = [x => 105, y => 105])
    state = TearingState(sys)
    @test_logs (:warn, r"state_priority") match_mode = :any ModelingToolkit.eliminate_perfect_aliases!(state)
    # Exactly one of the two variables is eliminated as observed.
    n_elim = count(
        eq -> isequal(eq.lhs, unwrap(x)) || isequal(eq.lhs, unwrap(y)),
        state.additional_observed
    )
    @test n_elim == 1
end

@testset "Perfect aliases do not eliminate irreducible variables" begin
    @variables x(t) y(t)
    @variables e(t) [irreducible = true]
    @variables c(t) [irreducible = true] d(t) [irreducible = true]
    # Two independent alias groups:
    #   * {x, e}     -- one irreducible; the non-irreducible `x` is eliminated as observed
    #   * {c, d, y}  -- two irreducibles + one non-irreducible. `y` is eliminated, both
    #                   irreducibles remain unknowns, bound by the surviving alias
    #                   equation between them.
    @mtkcompile sys = System(
        [
            D(x) ~ x,
            D(c) ~ -c,
            e ~ x,
            c ~ d,
            y ~ c,
        ], t
    )

    @test Set(unknowns(sys)) == Set([e, c, d])
end

@testset "`eliminate_perfect_aliases!` correctly handles unscalarized arrays" begin
    @variables x(t)[1:2] y(t)[1:2]
    @mtkcompile sys = System([D(x) ~ x, y[1] ~ y[2], dot(x, y) ~ 1], t)
    @test any(equations(sys)) do eq
        isequal(eq, 0 ~ 1 - dot(x, Symbolics.SConst([y[2], y[2]]))) ||
            isequal(eq, 0 ~ 1 - dot(x, Symbolics.SConst([y[1], y[1]])))
    end
end

@testset "Perfect aliases detect negated form `x ~ -y`" begin
    # Run only `eliminate_perfect_aliases!` on a fresh `TearingState`, so we
    # observe the output of *this* pass rather than the cumulative effect of
    # the rest of `mtkcompile`. The pass mutates `state` in place:
    #   - `state.additional_observed` collects `v ~ rhs` for each eliminated v
    #   - surviving equations of the system are returned by `equations(state.sys)`
    #   - irreducibles forced to zero by sign conflicts have one alias equation
    #     rewritten in place to `irr ~ 0`.

    # Numerically evaluate the RHS of an observed equation at a chosen value
    # of the surviving unknown, without depending on whether the symbolic
    # representation of `-x` is `-x` or `-1*x` etc.
    eval_rhs(rhs, var, val) = Symbolics.value(Symbolics.substitute(rhs, Dict(var => val)))

    # Find the entry of `state.additional_observed` whose lhs matches `v`.
    obs_for(state, v) = only(filter(eq -> isequal(eq.lhs, v), state.additional_observed))

    @testset "basic negated alias `y ~ -x`" begin
        @variables x(t) y(t)
        @named sys = System([D(x) ~ -x, y ~ -x], t; state_priorities = [x => 10])
        state = TearingState(sys)
        ModelingToolkit.eliminate_perfect_aliases!(state)
        # y ~ -x recorded as observed
        @test eval_rhs(obs_for(state, y).rhs, x, 3.0) == -3.0
        # The alias equation is gone; the dynamics `D(x) ~ -x` is the only one left.
        eqs = equations(state.sys)
        @test length(eqs) == 1
        @test isequal(only(eqs).lhs, D(x))
        # `y` is no longer in `fullvars` after `rm_eqs_vars!`.
        @test !any(isequal(unwrap(y)), state.fullvars)
    end

    @testset "mixed chain `x ~ y, y ~ -z`" begin
        @variables x(t) y(t) z(t)
        @named sys = System(
            [D(x) ~ -x, x ~ y, y ~ -z], t; state_priorities = [x => 10]
        )
        state = TearingState(sys)
        ModelingToolkit.eliminate_perfect_aliases!(state)
        # y ~ x (sign +1) and z ~ -x (sign -1: y has +1, z has -1 via `y ~ -z`)
        @test eval_rhs(obs_for(state, y).rhs, x, 3.0) == 3.0
        @test eval_rhs(obs_for(state, z).rhs, x, 3.0) == -3.0
    end

    @testset "double-negation transitivity `x ~ -y, y ~ -z` ⇒ `x ~ z`" begin
        @variables x(t) y(t) z(t)
        @named sys = System(
            [D(x) ~ -x, x ~ -y, y ~ -z], t; state_priorities = [x => 10]
        )
        state = TearingState(sys)
        ModelingToolkit.eliminate_perfect_aliases!(state)
        # y ~ -x (sign -1), z ~ x (sign +1: two negations cancel)
        @test eval_rhs(obs_for(state, y).rhs, x, 3.0) == -3.0
        @test eval_rhs(obs_for(state, z).rhs, x, 3.0) == 3.0
    end

    @testset "negated alias with irreducible target" begin
        @variables x(t) [irreducible = true]
        @variables y(t)
        @named sys = System([D(x) ~ -x, y ~ -x], t)
        state = TearingState(sys)
        ModelingToolkit.eliminate_perfect_aliases!(state)
        # x is irreducible, so y is the eliminated one: y ~ -x is observed.
        @test eval_rhs(obs_for(state, y).rhs, x, 2.5) == -2.5
        # x must be marked `always_present` by the pass (so downstream keeps it).
        x_idx = findfirst(isequal(unwrap(x)), state.fullvars)
        @test x_idx !== nothing && state.always_present[x_idx]
    end

    @testset "conflicting signs force both to zero (no irreducibles)" begin
        @variables x(t) y(t) z(t)
        # `x ~ y` and `x ~ -y` together imply `x = y = 0`. Neither is irreducible,
        # so both are substituted to `0` (no alias-equation pin needed) and both
        # alias equations become `0 ~ 0` and are removed.
        @named sys = System([D(z) ~ z, x ~ y, x ~ -y], t)
        state = TearingState(sys)
        ModelingToolkit.eliminate_perfect_aliases!(state)
        @test iszero(Symbolics.value(obs_for(state, x).rhs))
        @test iszero(Symbolics.value(obs_for(state, y).rhs))
        # Both original alias equations are removed; only `D(z) ~ z` survives.
        eqs = equations(state.sys)
        @test length(eqs) == 1
        @test isequal(only(eqs).lhs, D(z))
    end

    @testset "conflict with 3 irreducibles, one only on a single eq" begin
        @variables a(t) [irreducible = true]
        @variables b(t) [irreducible = true]
        @variables c(t) [irreducible = true]
        @variables w(t)
        # All three are irreducible, all in one conflict group via `a ~ b`,
        # `a ~ c`, `a ~ -c`. Each irreducible needs its own `irr ~ 0` equation.
        # `b` only appears on `e1 = a ~ b`. A previous greedy claiming
        # algorithm (prefer v1, fall back to v2) would have stranded `b`:
        # e1 → claim a, e2=(a,c) → claim c, e3=(a,c) → both claimed, drop. So
        # b would survive as an unconstrained unknown. With the matching-free
        # rewrite we now do (any irr ↔ any group eq), all three irreducibles
        # get their own pinned equation.
        @named sys = System([D(w) ~ w, a ~ b, a ~ c, a ~ -c], t)
        state = TearingState(sys)
        ModelingToolkit.eliminate_perfect_aliases!(state)
        eqs = equations(state.sys)
        for v in (a, b, c)
            @test any(
                eq -> isequal(eq.lhs, unwrap(v)) && iszero(Symbolics.value(eq.rhs)), eqs
            )
        end
        @test any(eq -> isequal(eq.lhs, D(w)), eqs)
    end

    @testset "conflicting signs with irreducible: alias eq rewritten to `v ~ 0`" begin
        @variables x(t) [irreducible = true]
        @variables y(t) z(t)
        # `x ~ y` and `x ~ -y` with `x` irreducible: `x` cannot be removed, so one
        # alias equation is rewritten to `x ~ 0` and the other is removed.
        @named sys = System([D(z) ~ z, x ~ y, x ~ -y], t)
        state = TearingState(sys)
        ModelingToolkit.eliminate_perfect_aliases!(state)
        # y is non-irreducible: substituted to 0 in observed.
        @test iszero(Symbolics.value(obs_for(state, y).rhs))
        # x is NOT in additional_observed (it survives as an unknown).
        @test !any(eq -> isequal(eq.lhs, unwrap(x)), state.additional_observed)
        # The system carries `x ~ 0` as a real equation now, alongside the
        # untouched `D(z) ~ z`.
        eqs = equations(state.sys)
        @test length(eqs) == 2
        @test any(eq -> isequal(eq.lhs, unwrap(x)) && iszero(Symbolics.value(eq.rhs)), eqs)
        @test any(eq -> isequal(eq.lhs, D(z)), eqs)
    end
end

@testset "Substitution rebuilds equation incidence" begin
    # Calls `find_perfect_aliases!` directly so we can inspect the bipartite
    # graph *before* `rm_eqs_vars!` renumbers it.
    sneighbors = ModelingToolkit.BipartiteGraphs.𝑠neighbors

    @testset "zero substitution drops multiplicatively annihilated var" begin
        @variables x(t) y(t) a(t) b(t)
        # `x ~ y` and `x ~ -y` is a sign conflict ⇒ x = y = 0. The non-alias
        # eq `D(a) ~ x*b + a` simplifies to `D(a) ~ a` after `x → 0`, so `b`
        # is annihilated and its edge to that eq must be dropped even though
        # `b` itself wasn't substituted.
        @named sys = System(
            [D(a) ~ x * b + a, D(b) ~ b, x ~ y, x ~ -y], t
        )
        state = TearingState(sys)

        eqs_before = collect(ModelingToolkit.equations(state))
        da_ieq = findfirst(eq -> isequal(eq.lhs, D(a)), eqs_before)
        b_idx = findfirst(isequal(unwrap(b)), state.fullvars)
        # Precondition: b is incident on the D(a) eq before the pass runs.
        @test b_idx in sneighbors(state.structure.graph, da_ieq)

        ModelingToolkit.find_perfect_aliases!(state, Int[], Int[])

        @test !(b_idx in sneighbors(state.structure.graph, da_ieq))
    end

    @testset "alias substitution drops cancelled target var" begin
        @variables x(t) y(t) w(t)
        # `x ~ y` is a consistent alias. `x` gets eliminated in favor of `y`.
        # The non-alias eq `D(w) ~ x - y + w` becomes `D(w) ~ w`, so `y` is no
        # longer in the equation even though it was the alias *target*.
        @named sys = System(
            [D(w) ~ x - y + w, x ~ y], t;
            state_priorities = [y => 10]
        )
        state = TearingState(sys)

        eqs_before = collect(ModelingToolkit.equations(state))
        dw_ieq = findfirst(eq -> isequal(eq.lhs, D(w)), eqs_before)
        y_idx = findfirst(isequal(unwrap(y)), state.fullvars)
        @test y_idx in sneighbors(state.structure.graph, dw_ieq)

        ModelingToolkit.find_perfect_aliases!(state, Int[], Int[])

        @test !(y_idx in sneighbors(state.structure.graph, dw_ieq))
    end
end

@testset "`eliminate_perfect_aliases!` correctly handles unscalarized array variables" begin
    @variables x(t)[1:2] y(t) z(t) [state_priority = -10] w(t)
    @named sys = System([0 ~ dot(x, Symbolics.SConst([y, z])), z ~ w], t)
    ts = TearingState(sys)
    ModelingToolkit.eliminate_perfect_aliases!(ts)
    # `z ~ w` is eliminated, `w` is substituted into the remaining equation
    @test length(equations(ts)) == 1
    @test length(ts.fullvars) == 4
    # Prior to the fix, this only had `y` and `w` since unscalarized array incidence
    # was not handled correctly in the pass.
    @test issetequal(𝑠neighbors(ts.structure.graph, 1), 1:4)
end

@testset "`eliminate_zero_variables!`" begin
    SS = ModelingToolkit.StateSelection
    analytically_integrated = ModelingToolkit.analytically_integrated

    # Build a fresh `TearingState`, populate the integer-coefficient linear
    # incidence matrix `mm` via `linear_subsys_adjmat!`, and run the zero-variable
    # elimination pass directly on it. Mirrors how `mtkcompile!` invokes the pass.
    function run_zero_var_pass(sys; kwargs...)
        ts = TearingState(sys)
        SS.complete!(ts.structure)
        mm = SS.linear_subsys_adjmat!(ts)
        mm, modified = ModelingToolkit.eliminate_zero_variables!(ts, mm; kwargs...)
        return ts, mm, modified
    end

    # Is variable `v` still present in `fullvars`?
    hasvar(ts, v) = any(isequal(unwrap(v)), ts.fullvars)
    # The `additional_observed` equation whose lhs matches the (already unwrapped) `v`.
    obs_for(ts, v) = only(filter(o -> isequal(o.lhs, v), ts.additional_observed))
    # `default_toterm` of a differential variable, e.g. `D(x) -> xˍt(t)`.
    tt(v) = ModelingToolkit.default_toterm(unwrap(v))
    obs_lhss(ts) = [o.lhs for o in ts.additional_observed]
    all_obs_zero(ts) = all(o -> iszero(Symbolics.value(o.rhs)), ts.additional_observed)
    # Parameters the pass introduced (the `t = 0` values of analytically integrated
    # variables), i.e. those absent from the original system `sys`.
    new_params(ts, sys) = setdiff(parameters(ts.sys), parameters(sys))

    @testset "a purely algebraic zero variable is eliminated" begin
        # `0 ~ 2y` forces `y == 0`. The row for that equation in `mm` has a single
        # non-zero entry, which is how the pass identifies `y` as zero.
        @variables x(t) y(t)
        @named sys = System([D(x) ~ x + y, 0 ~ 2y], t)
        ts, mm, modified = run_zero_var_pass(sys)
        @test modified
        @test !hasvar(ts, y)
        @test hasvar(ts, x) && hasvar(ts, D(x))
        # `y ~ 0` is recorded as an observed equation.
        @test isequal(only(ts.additional_observed).lhs, unwrap(y))
        @test iszero(Symbolics.value(only(ts.additional_observed).rhs))
        # `y` is not a derivative, so nothing is analytically integrated and no
        # `t = 0` parameters are introduced.
        @test isempty(analytically_integrated(ts.sys))
        @test isempty(new_params(ts, sys))
    end

    @testset "returns `mm` unchanged when there are no zero variables" begin
        @variables x(t) y(t)
        @named sys = System([D(x) ~ y, D(y) ~ -x], t)
        ts, mm, modified = run_zero_var_pass(sys)
        @test !modified
        @test isempty(ts.additional_observed)
        # `D(x), x, D(y), y` all survive.
        @test length(ts.fullvars) == 4
    end

    @testset "derivatives of a zero variable are also zeroed" begin
        # `0 ~ 2x` forces `x == 0`, and hence `D(x) == 0`. Both must be removed and
        # recorded as observed to be zero.
        @variables x(t)
        @named sys = System([0 ~ 2x, D(x) ~ 3x], t, [x], [])
        ts, mm, modified = run_zero_var_pass(sys)
        @test modified
        @test !hasvar(ts, x) && !hasvar(ts, D(x))
        @test issetequal(obs_lhss(ts), [unwrap(x), tt(D(x))])
        @test all_obs_zero(ts)
        # Every variable here is genuinely zero (not integrated), so nothing is
        # analytically integrated and no `t = 0` parameters are introduced.
        @test isempty(analytically_integrated(ts.sys))
        @test isempty(new_params(ts, sys))
    end

    @testset "a zero first derivative integrates to a constant" begin
        # `D(x) ~ 0` means `x` is constant. Rather than remaining an unknown, `x` is
        # "analytically integrated": eliminated in favour of `x ~ x0`, where `x0` is a
        # fresh parameter holding its value at `t = 0`.
        @variables x(t)
        @named sys = System([D(x) ~ 0], t, [x], [])
        ts, mm, modified = run_zero_var_pass(sys)
        @test modified
        # `D(x) ~ 0`
        @test iszero(Symbolics.value(obs_for(ts, tt(D(x))).rhs))
        # `x` is recorded as analytically integrated.
        @test issetequal(keys(analytically_integrated(ts.sys)), [unwrap(x)])
        # Exactly one fresh `t = 0` parameter is introduced, and `x ~ x0`.
        x0 = only(new_params(ts, sys))
        @test isequal(obs_for(ts, unwrap(x)).rhs, x0)
        # The `t = 0` parameter is solvable (bound to `missing`).
        @test isequal(ModelingToolkit.get_bindings(ts.sys)[x0], ModelingToolkit.COMMON_MISSING)
        # The analytical derivative of `x` is recorded as zero.
        @test iszero(Symbolics.value(ts.analytical_derivatives[unwrap(D(x))]))
    end

    @testset "a zero second derivative integrates to a polynomial in time" begin
        # `D(D(x)) ~ 0` is constant-velocity motion. Neither `x` nor `D(x)` is zero;
        # both are analytically integrated into polynomials in `t`, parameterised by
        # their `t = 0` values `x0` and `v0`:
        #   D(x) ~ v0
        #   x    ~ x0 + v0 * t
        @variables x(t)
        @named sys = System([D(D(x)) ~ 0], t, [x], [])
        ts, mm, modified = run_zero_var_pass(sys)
        @test modified
        # `D(D(x)) ~ 0`
        @test iszero(Symbolics.value(obs_for(ts, tt(D(D(x)))).rhs))
        # Both `x` and `D(x)` are analytically integrated, introducing two fresh
        # `t = 0` parameters.
        @test issetequal(keys(analytically_integrated(ts.sys)), [unwrap(x), unwrap(D(x))])
        newps = new_params(ts, sys)
        @test length(newps) == 2
        # `D(x) ~ v0`, a single `t = 0` parameter.
        v0 = obs_for(ts, tt(D(x))).rhs
        @test any(isequal(v0), newps)
        # `x ~ x0 + v0 * t`: the value at `t = 0` and the coefficient of `t` are exactly
        # the two new parameters, and the linear coefficient is `v0`.
        rhs_x = obs_for(ts, unwrap(x)).rhs
        x0 = Symbolics.value(Symbolics.substitute(rhs_x, Dict(unwrap(t) => 0)))
        dxdt = (rhs_x - x0) / t
        @test issetequal([x0, dxdt], newps)
        @test isequal(dxdt, v0)
        # Evaluate the polynomial at `x0 = 2`, `v0 = 3`, `t = 5`: `2 + 3*5 = 17`.
        val = Symbolics.value(Symbolics.substitute(rhs_x, Dict(x0 => 2.0, v0 => 3.0, unwrap(t) => 5.0)))
        @test val == 17.0
        # Analytical derivatives are recorded for reassembly.
        @test iszero(Symbolics.value(ts.analytical_derivatives[D(tt(D(x)))]))
        @test isequal(ts.analytical_derivatives[unwrap(D(x))], v0)
    end

    @testset "cascaded zeros are found within a single pass" begin
        # `0 ~ 2z` gives `z == 0`. Substituting into `0 ~ y - z` leaves `0 ~ y`, so
        # `y == 0` too. Both are discovered in one invocation because the pass queues
        # equations incident on a freshly-zeroed variable.
        @variables x(t) y(t) z(t)
        @named sys = System([0 ~ 2z, 0 ~ y - z, D(x) ~ x + y + z], t, [x, y, z], [])
        ts, mm, modified = run_zero_var_pass(sys)
        @test modified
        @test !hasvar(ts, y) && !hasvar(ts, z)
        @test hasvar(ts, x)
        @test issetequal(obs_lhss(ts), [unwrap(y), unwrap(z)])
        @test all_obs_zero(ts)
    end

    @testset "`eliminate_zero_variables_fixpoint!` discovers zeros across iterations" begin
        # `w` is a parameter, so `0 ~ w*z + y` is not integer-linear and is absent from
        # `mm` initially. Only after `z -> 0` (first iteration) does it become `0 ~ y`
        # and get added to `mm`, so a second iteration is required to find `y == 0`.
        @variables x(t) y(t) z(t)
        @parameters w
        @named sys = System([0 ~ 2z, 0 ~ w * z + y, D(x) ~ x + y + z], t, [x, y, z], [w])

        # A single pass only finds `z`.
        ts1 = TearingState(sys)
        SS.complete!(ts1.structure)
        mm1 = SS.linear_subsys_adjmat!(ts1)
        mm1, _ = ModelingToolkit.eliminate_zero_variables!(ts1, mm1)
        @test !hasvar(ts1, z)
        @test hasvar(ts1, y)

        # The fixpoint driver iterates until convergence and finds both.
        ts2 = TearingState(sys)
        SS.complete!(ts2.structure)
        mm2 = SS.linear_subsys_adjmat!(ts2)
        mm2 = ModelingToolkit.eliminate_zero_variables_fixpoint!(ts2, mm2)
        @test !hasvar(ts2, y) && !hasvar(ts2, z)
        @test issetequal(obs_lhss(ts2), [unwrap(y), unwrap(z)])

        # `maxiters = 1` stops after the first iteration, matching the single pass.
        ts3 = TearingState(sys)
        SS.complete!(ts3.structure)
        mm3 = SS.linear_subsys_adjmat!(ts3)
        mm3 = ModelingToolkit.eliminate_zero_variables_fixpoint!(ts3, mm3; maxiters = 1)
        @test hasvar(ts3, y)
    end

    @testset "zero variables are substituted through `mtkcompile`" begin
        # End-to-end: `0 ~ 2y` should vanish and `y -> 0` should be substituted into
        # the surviving dynamics, leaving `D(x) ~ x`.
        @variables x(t) y(t)
        @mtkcompile sys = System([D(x) ~ x + y, 0 ~ 2y], t)
        @test length(equations(sys)) == 1
        eq = only(equations(sys))
        @test isequal(eq.lhs, unwrap(D(x)))
        @test isequal(eq.rhs, unwrap(x))
        @test any(o -> isequal(o.lhs, unwrap(y)) && iszero(Symbolics.value(o.rhs)), observed(sys))
    end

    @testset "analytically integrated variables reconstruct the trajectory" begin
        # `D(D(x)) ~ 0` compiles to a system with no equations and no unknowns: the
        # entire trajectory lives in `observed` as `x ~ x0 + v0 * t`, with `x0` and `v0`
        # parameters standing in for the `t = 0` position and velocity.
        @variables x(t)
        @mtkcompile sys = System([D(D(x)) ~ 0], t, [x], [])
        @test isempty(equations(sys))
        @test isempty(unknowns(sys))
        @test issetequal(keys(analytically_integrated(sys)), [unwrap(x), unwrap(D(x))])

        # The observed expression for `x` is `x0 + v0 * t`, and its two parameters are
        # exactly the parameters of the compiled system.
        obsx = only(filter(o -> isequal(o.lhs, unwrap(x)), observed(sys)))
        x0 = Symbolics.value(Symbolics.substitute(obsx.rhs, Dict(unwrap(t) => 0)))
        v0 = (obsx.rhs - x0) / t
        @test issetequal([x0, v0], parameters(sys))

        # Initial conditions for `x` and `D(x)` determine `x0` and `v0`; the observed
        # function then reconstructs the trajectory (`x(0) = 2`, `x(5) = 2 + 3*5 = 17`).
        prob = ODEProblem(sys, [x => 2.0, D(x) => 3.0], (0.0, 5.0))
        @test prob.f.observed(unwrap(x), prob.u0, prob.p, 0.0) == 2.0
        @test prob.f.observed(unwrap(x), prob.u0, prob.p, 5.0) == 17.0
        # With no initialization information at all, construction fails.
        @test_throws ModelingToolkit.IncompleteInitializationError ODEProblem(sys, [], (0.0, 5.0))
    end
end
