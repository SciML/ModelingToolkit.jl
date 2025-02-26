using ModelingToolkit, Test
using ModelingToolkit: t_nounits as t
using StableRNGs

k = ShiftIndex(t)
rng = StableRNG(22525) 

# Shift(t, -1)(x(t)) - x_{t-1}(t)
# -3 - x(t) + x(t)*x_{t-1}
@testset "Correct ImplicitDiscreteFunction" begin
    @variables x(t) = 1
    @mtkbuild sys = ImplicitDiscreteSystem([x(k) ~ x(k)*x(k-1) - 3], t)
    tspan = (0, 10)
    f = ImplicitDiscreteFunction(sys)
    u_next = [3., 1.5]
    @test f(u_next, [2.,3.], [], t) ≈ [0., 0.]
    u_next = [0., 0.]
    @test f(u_next, [2.,3.], [], t) ≈ [3., -3.]
    
    resid = rand(2)
    f(resid, u_next, [2.,3.], [], t)
    @test resid ≈ [3., -3.]
    
    prob = ImplicitDiscreteProblem(sys, [x(k-1) => 3.], tspan)
    @test prob.u0 == [3., 1.]
    prob = ImplicitDiscreteProblem(sys, [], tspan)
    @test prob.u0 == [1., 1.]
    @variables x(t)
    @mtkbuild sys = ImplicitDiscreteSystem([x(k) ~ x(k)*x(k-1) - 3], t)
    @test_throws ErrorException prob = ImplicitDiscreteProblem(sys, [], tspan)
end

# Test solvers
@testset "System with algebraic equations" begin
    @variables x(t) y(t)
    eqs = [x(k) ~ x(k-1) + x(k-2), 
           x^2 ~ 1 - y^2]
    @mtkbuild sys = ImplicitDiscreteSystem(eqs, t)
    f = ImplicitDiscreteFunction(sys)

    function correct_f(u_next, u, p, t) 
        [u[2] - u_next[1],
         u[1] + u[2] - u_next[2],
         1 - (u_next[1]+u_next[2])^2 - u_next[3]^2]
    end

    for _ in 1:10
        u_next = rand(rng, 3)
        u = rand(rng, 3)
        @test correct_f(u_next, u, [], 0.) ≈ f(u_next, u, [], 0.)
    end

    # Initialization is satisfied.
    prob = ImplicitDiscreteProblem(sys, [x(k-1) => 3.], tspan)
    @test (prob.u0[1] + prob.u0[2])^2 + prob.u0[3]^2 ≈ 1
end

@testset "System with algebraic equations, implicit difference equations, explicit difference equations" begin
    @variables x(t) y(t)
    eqs = [x(k) ~ x(k-1) + x(k-2),
           y(k) ~ x(k) + x(k-2)*y(k-1)]
    @mtkbuild sys = ImplicitDiscreteSystem(eqs, t)
end
