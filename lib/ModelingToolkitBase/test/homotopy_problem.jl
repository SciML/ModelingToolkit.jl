using ModelingToolkitBase
using ModelingToolkitBase: homotopy, has_any_homotopy, get_nonlinear_problem_type,
    t_nounits as t, D_nounits as D
using SciMLBase
using NonlinearSolve   # HomotopySweep (re-exported) + inner nonlinear solvers
using Test

@testset "direct HomotopyProblem(sys, op): construction + residual endpoints" begin
    @variables y
    @mtkcompile sys = System([0 ~ homotopy(atan(y - 3), y)])
    prob = HomotopyProblem(sys, [y => 12.0])
    @test prob isa SciMLBase.HomotopyProblem
    @test prob.λspan == (0.0, 1.0)                        # inherited SciMLBase default
    @test !hasfield(typeof(prob), :homotopy_parameter)    # merged API: λ is an arg, not a field
    # the swept residual genuinely depends on λ (endpoints differ through the wrap)
    @test !(prob.f(copy(prob.u0), prob.p, 0.0) ≈ prob.f(copy(prob.u0), prob.p, 1.0))
    # a custom λspan threads through to the constructed problem
    @test HomotopyProblem(sys, [y => 12.0]; λspan = (0.0, 0.5)).λspan == (0.0, 0.5)
    # expression = Val{true} is an explicit (documented) error, not silent misgen
    @test_throws ArgumentError HomotopyProblem(sys, [y => 12.0]; expression = Val{true})

    u = copy(prob.u0)
    p = prob.p
    @test only(prob.f(u, p, 0.0)) ≈ 12.0 atol = 1.0e-10              # simplified endpoint
    @test only(prob.f(u, p, 1.0)) ≈ atan(9.0) atol = 1.0e-10        # actual endpoint
    @test only(prob.f(u, p, 0.5)) ≈ 0.5 * 12.0 + 0.5 * atan(9.0) atol = 1.0e-10  # convex blend
    # in-place residual variant f(du, u, p, λ): same wrapper, arity 4
    if SciMLBase.isinplace(prob)
        du = similar(u)
        prob.f(du, u, p, 0.5)
        @test only(du) ≈ only(prob.f(u, p, 0.5)) atol = 1.0e-10
    end
end

@testset "AbstractNonlinearProblem(sys, op) selects HomotopyProblem ⟺ homotopy present" begin
    @variables y
    @mtkcompile hsys = System([0 ~ homotopy(atan(y - 3), y)])
    @test get_nonlinear_problem_type(hsys) === SciMLBase.HomotopyProblem
    @test AbstractNonlinearProblem(hsys, [y => 12.0]) isa SciMLBase.HomotopyProblem

    @variables z
    @mtkcompile nsys = System([0 ~ atan(z - 3)])
    @test get_nonlinear_problem_type(nsys) === SciMLBase.NonlinearProblem
    np = AbstractNonlinearProblem(nsys, [z => 12.0])
    @test np isa SciMLBase.NonlinearProblem
    @test !(np isa SciMLBase.HomotopyProblem)

    # HomotopyProblem on a system without `homotopy` is a loud error, not a
    # silent λ-independent problem.
    @test_throws ArgumentError HomotopyProblem(nsys, [z => 12.0])
end

@testset "constructor converts a time-dependent system (ODE -> NonlinearSystem)" begin
    # A time-dependent system exercises the is_time_dependent -> NonlinearSystem
    # conversion in the constructor (D(x) substituted to 0), which must still leave
    # the homotopy node intact and produce a swept HomotopyProblem.
    @variables x(t) y(t)
    @mtkcompile sys = System(
        [D(x) ~ -x, 0 ~ homotopy(atan(y - 3), y)], t; guesses = [y => 12.0]
    )
    prob = HomotopyProblem(sys, [x => 1.0, y => 12.0])
    @test prob isa SciMLBase.HomotopyProblem
    # residual is λ-dependent through the conversion + lowering
    u = copy(prob.u0)
    @test !(prob.f(u, prob.p, 0.0) ≈ prob.f(u, prob.p, 1.0))
    sol = solve(prob)
    @test SciMLBase.successful_retcode(sol)
    @test sol[y] ≈ 3.0 atol = 1.0e-6
end

@testset "e2e rescue: continuation walks an out-of-basin guess to the true root" begin
    @variables y
    @mtkcompile sys = System([0 ~ homotopy(atan(y - 3), y)])
    prob = HomotopyProblem(sys, [y => 12.0])
    # a `HomotopyProblem` with no algorithm defaults to the continuation sweep
    sol = solve(prob)
    @test SciMLBase.successful_retcode(sol)
    @test sol[y] ≈ 3.0 atol = 1.0e-6
    # an explicit `HomotopySweep` gives the same answer
    sol2 = solve(prob, HomotopySweep())
    @test sol2[y] ≈ 3.0 atol = 1.0e-6
end

@testset "structured p passes through the sweep untouched" begin
    @variables y
    @parameters c = 3.0
    @mtkcompile sys = System([0 ~ homotopy(atan(y - c), y)])
    prob = HomotopyProblem(sys, [y => 12.0, c => 3.0])
    @test prob.p isa ModelingToolkitBase.MTKParameters
    sol = solve(prob)
    @test sol[y] ≈ 3.0 atol = 1.0e-6
    @test sol.ps[c] == 3.0                # user parameter unscathed by the sweep
end

@testset "two homotopy() calls share the single sweep" begin
    @variables y z
    @mtkcompile sys = System(
        [0 ~ homotopy(atan(y - 3), y), 0 ~ homotopy(atan(z + 2), z)]
    )
    prob = HomotopyProblem(sys, [y => 12.0, z => -12.0])
    @test prob isa SciMLBase.HomotopyProblem
    sol = solve(prob)
    @test sol[y] ≈ 3.0 atol = 1.0e-6
    @test sol[z] ≈ -2.0 atol = 1.0e-6
end
