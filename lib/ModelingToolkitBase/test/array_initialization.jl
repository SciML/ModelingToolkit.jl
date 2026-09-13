using ModelingToolkitBase, SciMLBase, Test

@testset "Array derivative initialization guesses" begin
    @parameters t
    D = Differential(t)
    for dims in ((3,), (2, 3))
        @testset "shape $dims" begin
            ax = map(Base.OneTo, dims)
            @variables x(t)[ax...]
            values = reshape(collect(1.0:prod(dims)), dims)
            @named raw = System([D(x) ~ -x], t, vec(collect(x)), [])
            sys = complete(raw)
            op = [x => values, D(x) => -values]
            prob = DAEProblem(sys, op, (0.0, 1.0))
            @test prob.u0 == vec(values)
            @test prob.du0 == -vec(values)
            residual = similar(prob.u0)
            prob.f(residual, prob.du0, prob.u0, prob.p, 0.0)
            @test residual == zeros(length(values))
        end
    end
end
