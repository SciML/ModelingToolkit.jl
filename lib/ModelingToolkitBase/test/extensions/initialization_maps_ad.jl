using ForwardDiff, ModelingToolkitBase, Test, Zygote
using ModelingToolkitBase: t_nounits as t, D_nounits as D, SciMLBase
using SymbolicIndexingInterface: ProblemState

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

@testset "FullSpecialize state map gradients" begin
    @variables x(t) y(t)
    @mtkcompile sys = System(
        [D(x) ~ -x, D(y) ~ -y], t;
        initialization_eqs = [x^3 + x ~ 2, y ~ 2x + 1]
    )
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(
        sys, [], (0.0, 1.0); guesses = [x => 1.0, y => 1.0]
    )
    data = prob.f.initialization_data
    initprob = data.initializeprob
    loss(u) = sum(abs2, data.initializeprobmap(ProblemState(; u, p = initprob.p, t = 0.0)))
    @test only(Zygote.gradient(loss, initprob.u0)) ≈ ForwardDiff.gradient(loss, initprob.u0)
end
