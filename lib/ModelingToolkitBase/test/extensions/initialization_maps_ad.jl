using ForwardDiff, ModelingToolkitBase, Test, Zygote
using ModelingToolkitBase: t_nounits as t, D_nounits as D, SciMLBase
using SymbolicIndexingInterface: ProblemState
using SciMLStructures: Tunable, replace, canonicalize
using SymbolicIndexingInterface: parameter_values, parameter_index, remake_buffer
using ChainRulesCore: NoTangent
import ChainRulesCore

@testset "FullSpecialize parameter buffer gradients" begin
    for prototype in ([0.0], ModelingToolkitBase.SVector(0.0), Any[0.0])
        loss(x) = sum(abs2, ModelingToolkitBase._static_initialization_buffer(prototype, (x,)))
        @test only(Zygote.gradient(loss, 2.0)) == 4.0
    end

    function tagged_loss(x)
        buffer = ModelingToolkitBase._static_initialization_buffer([Val(1)], (Val(1),))
        return x^2 + length(buffer)
    end
    @test only(Zygote.gradient(tagged_loss, 2.0)) == 4.0
end

@testset "FullSpecialize initialization map gradients" begin
    @variables x(t) y(t)
    @parameters a = 1.0
    @mtkcompile sys = System(
        [D(x) ~ -a * x, D(y) ~ -y], t;
        initialization_eqs = [x^3 + x ~ 2, y ~ 2x + 1]
    )
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(
        sys, [], (0.0, 1.0); guesses = [x => 1.0, y => 1.0]
    )
    data = prob.f.initialization_data
    initprob = data.initializeprob
    loss(u) = sum(abs2, data.initializeprobmap(ProblemState(; u, p = initprob.p, t = 0.0)))
    @test only(Zygote.gradient(loss, initprob.u0)) ≈ ForwardDiff.gradient(loss, initprob.u0)

    function parameter_loss(p)
        state = ProblemState(; u = initprob.u0, p = replace(Tunable(), initprob.p, p), t = 0.0)
        return sum(abs2, data.initializeprobpmap(prob, state).tunable)
    end
    p = initprob.p.tunable
    expected = ForwardDiff.gradient(parameter_loss, p)
    @test expected == 2p
    @test only(Zygote.gradient(parameter_loss, p)) ≈ expected
end

@testset "`nonnumeric` cotangent is `NoTangent` (matches the Enzyme inactive rule)" begin
    @parameters m = 1.5 name::String = "a"
    @variables s(t) = 1.0 v(t) = 1.0
    # `name` is unused in the equations, so list it explicitly to keep it
    sys = mtkcompile(System([D(s) ~ v, D(v) ~ (1 - v) / m], t, [s, v], [m, name]; name = :model))
    prob = ODEProblem(sys, [], (0.0, 1.0))
    p0 = parameter_values(prob)
    @test !isempty(p0.nonnumeric)

    # constructor rrule: a rebuild via `SciMLStructures.replace` must not produce a
    # tangent for the nonnumeric portion
    buf = canonicalize(Tunable(), p0)[1]
    g = Zygote.gradient(
        p -> sum(replace(Tunable(), p, buf).initials), p0
    )[1]
    @test g.initials == ones(length(p0.initials))
    @test g.nonnumeric isa NoTangent || g.nonnumeric === nothing

    # `remake_buffer` rrule
    idxs = [parameter_index(sys, m)]
    _, back = ChainRulesCore.rrule(remake_buffer, sys, p0, idxs, [2.0])
    dp = back(ChainRulesCore.Tangent{typeof(p0)}(; tunable = ones(length(p0.tunable))))[3]
    @test dp.nonnumeric isa NoTangent

    # the passthrough itself is non-differentiable
    w = getfield(p0, :nonnumeric)
    @test Zygote.gradient(w -> length(w.buffer), w)[1] === nothing
end
