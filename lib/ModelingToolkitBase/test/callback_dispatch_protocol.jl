using ModelingToolkitBase, SciMLBase, OrdinaryDiffEqDefault, Test
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using ModelingToolkitBase: VectorAffect, VectorOptionalAffect, EMPTY_AFFECT

# An integrator stand-in that logs which affects ran and whether the step was marked modified.
mutable struct LogIntegrator
    calls::Vector{Int}
    modified::Bool
end
LogIntegrator() = LogIntegrator(Int[], true)
SciMLBase.derivative_discontinuity!(i::LogIntegrator, b::Bool) = (i.modified = b; nothing)

logger(n) = i -> push!(i.calls, n)

# Each case is run with the affects in a tuple and in a `Vector{Any}`: both representations
# route events the same way.
astuple(xs) = Tuple(xs)
asvector(xs) = Any[xs...]

@testset "VectorAffect routing ($(nameof(container)))" for container in (astuple, asvector)
    affects = container([logger(1), nothing, EMPTY_AFFECT, logger(4)])
    negs = container([nothing, EMPTY_AFFECT, nothing, logger(40)])
    va = VectorAffect([1, 2, 3, 4], affects, negs)

    # positive edge of 1 and negative edge of 4 run; the `EMPTY_AFFECT` entries are skipped
    i = LogIntegrator()
    va(i, [1, -1, 1, -1])
    @test i.calls == [1, 40]
    @test i.modified

    # only `nothing`/`EMPTY_AFFECT` entries triggered: nothing runs and the step is not modified
    i = LogIntegrator()
    va(i, [-1, 1, -1, 0])
    @test isempty(i.calls)
    @test !i.modified

    # no edge at all
    i = LogIntegrator()
    va(i, [0, 0, 0, 0])
    @test isempty(i.calls)
    @test !i.modified

    # two equations of one callback (shared affect index) and a second callback, all at once:
    # affects run in equation order, each edge choosing its own direction's affect
    va2 = VectorAffect(
        [1, 1, 2], container([logger(1), logger(2)]), container([logger(10), logger(20)])
    )
    i = LogIntegrator()
    va2(i, [1, -1, -1])
    @test i.calls == [1, 10, 20]
    @test i.modified
    i = LogIntegrator()
    va2(i, [0, 1, 1])
    @test i.calls == [1, 2]
end

@testset "apply_nth_affect out of range" begin
    @test_throws BoundsError ModelingToolkitBase.apply_nth_affect(
        (logger(1), logger(2)), 3, LogIntegrator()
    )
    @test_throws BoundsError ModelingToolkitBase.apply_nth_affect(
        Any[logger(1), logger(2)], 3, LogIntegrator()
    )
end

@testset "VectorOptionalAffect ($(nameof(container)))" for container in (astuple, asvector)
    i = LogIntegrator()
    VectorOptionalAffect(container([logger(1), nothing, logger(2)]))(nothing, nothing, nothing, i)
    @test i.calls == [1, 2]
    i = LogIntegrator()
    VectorOptionalAffect(container([nothing, nothing]))(nothing, nothing, nothing, i)
    @test isempty(i.calls)
end

# x = sin(t). On (0, 4), x crosses 0.5 upward at π/6 and downward at 5π/6, and crosses -0.5
# downward at 7π/6. Each crossing increments its own counter, and the initialize and finalize
# affects of the two callbacks set their own flags.
@testset "Vector callback edges, initialize and finalize in a solve" begin
    @variables x(t) = 0.0 v(t) = 1.0
    @discretes up1(t) = 0.0 dn1(t) = 0.0 up2(t) = 0.0 dn2(t) = 0.0 ini(t) = 0.0 fin(t) = 0.0
    sa(eqs, ps) = ModelingToolkitBase.SymbolicAffect(eqs; discrete_parameters = ps)
    cb1 = ModelingToolkitBase.SymbolicContinuousCallback(
        [x ~ 0.5], sa([up1 ~ Pre(up1) + 1], [up1]);
        affect_neg = sa([dn1 ~ Pre(dn1) + 1], [dn1]),
        initialize = sa([ini ~ 1.0], [ini]),
    )
    cb2 = ModelingToolkitBase.SymbolicContinuousCallback(
        [x ~ -0.5], sa([up2 ~ Pre(up2) + 1], [up2]);
        affect_neg = sa([dn2 ~ Pre(dn2) + 1], [dn2]),
        finalize = sa([fin ~ 1.0], [fin]),
    )
    @named osc = System(
        [D(x) ~ v, D(v) ~ -x], t, [x, v], [up1, dn1, up2, dn2, ini, fin];
        continuous_events = [cb1, cb2]
    )
    sys = mtkcompile(osc)
    prob = ODEProblem(sys, [], (0.0, 4.0))
    cb = prob.kwargs[:callback]
    cb = cb isa SciMLBase.CallbackSet ? only(cb.continuous_callbacks) : cb
    @test cb isa SciMLBase.VectorContinuousCallback
    integ = init(prob, DefaultODEAlgorithm(); abstol = 1.0e-10, reltol = 1.0e-10)
    @test integ.ps[ini] == 1.0
    @test integ.ps[fin] == 0.0
    solve!(integ)
    sol = integ.sol
    @test SciMLBase.successful_retcode(sol)
    @test integ.ps[up1] == 1.0
    @test integ.ps[dn1] == 1.0
    @test integ.ps[up2] == 0.0
    @test integ.ps[dn2] == 1.0
    @test integ.ps[fin] == 1.0
    @test sol(1.0; idxs = up1) == 1.0 && sol(0.4; idxs = up1) == 0.0
    @test sol(2.0; idxs = dn1) == 0.0 && sol(2.7; idxs = dn1) == 1.0
    @test sol(3.6; idxs = dn2) == 0.0 && sol(3.7; idxs = dn2) == 1.0
    @test sol(4.0; idxs = x) ≈ sin(4.0) atol = 1.0e-6
end
