using ModelingToolkit, Test
using ModelingToolkit: t_nounits as t

struct ArrayCallable5378 end
(::ArrayCallable5378)(x) = [x, 2x, 3x]
struct ScalarCallable5378 end
(::ScalarCallable5378)(x, i) = i * x

@testset "Callable parameter projection regression" begin
    @testset "Array-valued callable parameters in nonlinear problems" begin
        @variables y(t) r(t)[1:3]
        @parameters (f::ArrayCallable5378)(..)[1:3]::Real = ArrayCallable5378()
        @named sys = System([0 ~ y - 5.0, r ~ f(y)], t, [y, r], [f])
        csys = mtkcompile(sys; fully_determined = true)
        prob = NonlinearProblem(csys, [])
        residual = zeros(1)
        prob.f(residual, [5.0], prob.p)
        @test residual ≈ [0.0]
        @test length(parameters(csys)) == 1
        @test [prob.f.observed(r[i], [5.0], prob.p) for i in 1:3] ≈ [5.0, 10.0, 15.0]
    end

    @testset "Scalar callable parameter control" begin
        @variables y(t) r(t)[1:3]
        # Scalar calls already work and must remain supported.
        @parameters (g::ScalarCallable5378)(..) = ScalarCallable5378()
        @named scalar = System(
            [0 ~ y - 5.0; [r[i] ~ g(y, i) for i in 1:3]], t, [y, r], [g]
        )
        compiled_scalar = mtkcompile(scalar; fully_determined = true)
        scalar_prob = NonlinearProblem(compiled_scalar, [])
        residual = zeros(1)
        scalar_prob.f(residual, [5.0], scalar_prob.p)
        @test residual ≈ [0.0]
        @test length(parameters(compiled_scalar)) == 1
        @test [scalar_prob.f.observed(r[i], [5.0], scalar_prob.p) for i in 1:3] ≈ [5.0, 10.0, 15.0]
    end

    @testset "Genuine name collisions are still rejected" begin
        @variables duplicate_name(t)
        variable = duplicate_name
        @parameters duplicate_name
        parameter = duplicate_name
        @test_throws ArgumentError System([0 ~ variable - parameter], t, [variable], [parameter]; name = :conflict)
    end
end
