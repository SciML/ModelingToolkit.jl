using ModelingToolkitBase, Test
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using ModelingToolkitBase: MissingTspanError, get_tspan, has_tspan
using SciMLBase
import DiffEqNoiseProcess
using Setfield: @set

@testset "`tspan` is stored on and retrieved from the system" begin
    @variables x(t) = 1.0
    @named sys = System([D(x) ~ -x], t; tspan = (0.0, 5.0))
    @test has_tspan(sys)
    @test get_tspan(sys) == (0.0, 5.0)

    nosys = System([D(x) ~ -x], t; name = :sys)
    @test get_tspan(nosys) === nothing

    # `flatten`/`copy`/`complete`/`mtkcompile` must not lose it.
    @test get_tspan(flatten(sys)) == (0.0, 5.0)
    @test get_tspan(copy(sys)) == (0.0, 5.0)
    @test get_tspan(complete(sys)) == (0.0, 5.0)
    @test get_tspan(mtkcompile(sys)) == (0.0, 5.0)

    # `isapprox` compares systems field by field, so it has to notice the timespan.
    # The two systems here differ in nothing else.
    @test !isapprox(nosys, @set nosys.tspan = (0.0, 5.0))
end

# Each entry builds a system that is valid for the given problem type, along with the
# operating point to construct it with. `JumpProblem` is absent on purpose: see the comment
# at the top of `src/problems/jumpproblem.jl` for why it keeps a required `tspan`.
function tspan_problem_cases(tspan)
    cases = Pair{Any, Any}[]

    @variables x(t) = 1.0
    @parameters τ = 3.0
    push!(
        cases,
        ODEProblem => (mtkcompile(System([D(x) ~ -x / τ], t; tspan, name = :ode)), [])
    )

    @variables y(t) = 1.0
    push!(
        cases,
        DAEProblem => (
            mtkcompile(System([D(y) ~ -y], t; tspan, name = :dae)),
            [D(y) => -1.0],
        )
    )

    @variables z(t) = 1.0
    @brownians b
    push!(
        cases,
        SDEProblem => (mtkcompile(System([D(z) ~ -z + b], t; tspan, name = :sde)), [])
    )

    @variables w(..)
    push!(
        cases,
        DDEProblem => (
            mtkcompile(System([D(w(t)) ~ -w(t - 0.1)], t; tspan, name = :dde)),
            [w(t) => 1.0],
        )
    )

    @variables v(..)
    @brownians bs
    push!(
        cases,
        SDDEProblem => (
            mtkcompile(System([D(v(t)) ~ -v(t - 0.1) + bs], t; tspan, name = :sdde)),
            [v(t) => 1.0],
        )
    )

    @variables u(t) = 1.0
    push!(
        cases,
        BVProblem => (mtkcompile(System([D(u) ~ -u], t; tspan, name = :bvp)), [])
    )

    k = ShiftIndex(t)
    @variables d(t) = 1.0
    push!(
        cases,
        DiscreteProblem => (
            mtkcompile(System([d(k) ~ 0.5 * d(k - 1)], t; tspan, name = :disc)),
            [],
        )
    )
    @variables e(t) = 1.0
    push!(
        cases,
        ImplicitDiscreteProblem => (
            mtkcompile(System([e(k) ~ e(k) * e(k - 1) + 1], t; tspan, name = :idisc)),
            [],
        )
    )

    return cases
end

@testset "Problem constructors default to the system's `tspan`" begin
    for (T, (sys, op)) in tspan_problem_cases((0.0, 5.0))
        @testset "$T" begin
            # Omitted `tspan` falls back to the system's.
            @test T(sys, op).tspan == (0.0, 5.0)
            # An explicitly passed `tspan` always wins.
            @test T(sys, op, (0.0, 2.0)).tspan == (0.0, 2.0)
            # The fallback also reaches the `{iip}` and `{iip, specialize}` forms.
            @test T{true}(sys, op).tspan == (0.0, 5.0)
            @test T{true, SciMLBase.FullSpecialize}(sys, op).tspan == (0.0, 5.0)
        end
    end
end

@testset "Omitting `tspan` for a system without one errors" begin
    for (T, (sys, op)) in tspan_problem_cases(nothing)
        @testset "$T" begin
            @test_throws MissingTspanError T(sys, op)
            @test T(sys, op, (0.0, 2.0)).tspan == (0.0, 2.0)
        end
    end
end

@testset "`MissingTspanError` message names the system and the fix" begin
    err = MissingTspanError(:mysys)
    msg = sprint(showerror, err)
    @test occursin("mysys", msg)
    @test occursin("tspan", msg)
end
