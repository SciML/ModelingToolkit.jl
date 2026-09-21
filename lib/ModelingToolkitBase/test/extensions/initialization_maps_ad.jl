using ForwardDiff, ModelingToolkitBase, Test, Zygote
import Functors
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

    # Reading a nonnumeric value whose contents Zygote could attach a tangent to must not
    # produce a structured `(buffer = ...)` tangent: SciMLSensitivity reconciles parameter
    # cotangents with Functors and expects `nothing` there.
    @parameters tup::Tuple{Float64, Float64} = (1.0, 2.0)
    sys2 = mtkcompile(System([D(s) ~ v, D(v) ~ (1 - v) / m], t, [s, v], [m, tup]; name = :model2))
    p1 = parameter_values(ODEProblem(sys2, [], (0.0, 1.0)))
    g2 = Zygote.gradient(p -> sum(p.tunable) + p.nonnumeric[1][1][1], p1)[1]
    @test g2.tunable == ones(length(p1.tunable))
    @test g2.nonnumeric === nothing

    # SciMLSensitivity walks a parameter object together with its cotangent through
    # `Functors.fmap`, and also walks the parameter object alone to allocate buffers. The
    # wrapper must survive both: pairing with the `nothing` cotangent of the nonnumeric
    # slot, and having children so single-argument walks do not recurse on it forever.
    # (Functors cannot pair a non-empty tuple with `nothing` at all, so the paired walk is
    # only exercised with an empty nonnumeric portion, as it was before the wrapper.)
    zeroed = Functors.fmap(x -> x isa AbstractArray{<:Number} ? zero(x) : x, p1)
    @test zeroed isa typeof(p1)
    @test iszero(zeroed.tunable)
    @test zeroed.nonnumeric == p1.nonnumeric
    sys3 = mtkcompile(System([D(s) ~ v, D(v) ~ (1 - v) / m], t, [s, v], [m]; name = :model3))
    p3 = parameter_values(ODEProblem(sys3, [], (0.0, 1.0)))
    @test isempty(p3.nonnumeric)
    Δ = (; tunable = ones(length(p3.tunable)), initials = nothing, discrete = nothing,
        constant = nothing, nonnumeric = nothing, caches = nothing)
    out = Functors.fmap((y, x) -> x === nothing ? y : x, p3, Δ)
    @test out isa typeof(p3)
    @test out.tunable == ones(length(p3.tunable))
    @test out.nonnumeric === ()
end
