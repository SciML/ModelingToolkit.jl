using ModelingToolkitBase
using ModelingToolkitBase: homotopy, has_any_homotopy, get_nonlinear_problem_type,
    get_initialization_problem_type, InitializationProblem,
    t_nounits as t, D_nounits as D
using SciMLBase
using NonlinearSolve   # HomotopySweep / HomotopyPolyAlgorithm (re-exported) + inner solvers
using Test

# When the (square) initialization system carries Modelica `homotopy(actual, simplified)`
# nodes, `InitializationProblem` routes to a `SciMLBase.HomotopyProblem` (solved by
# continuation), selected with the same dispatch as `AbstractNonlinearProblem(isys, op)`.
# Systems without `homotopy` are unaffected (still `NonlinearProblem`).

@testset "square homotopy init system routes to HomotopyProblem" begin
    @variables x(t) y(t)
    @mtkcompile sys = System(
        [D(x) ~ -x, 0 ~ homotopy(atan(y - 3), y)], t; guesses = [y => 12.0]
    )

    # x is pinned by u0 ⇒ the init system for y is square (1 eq, 1 unknown) and carries
    # the homotopy node, so the selected init problem type is HomotopyProblem.
    iprob = InitializationProblem(sys, 0.0, [x => 1.0]; warn_initialize_determined = false)
    @test iprob isa SciMLBase.HomotopyProblem

    # solving the init problem walks the out-of-basin guess y=12 to the true root y=3
    sol = solve(iprob; abstol = 1.0e-12, reltol = 1.0e-12)
    @test SciMLBase.successful_retcode(sol)
    @test sol[y] ≈ 3.0 atol = 1.0e-6
end

@testset "non-homotopy init system is unaffected (NonlinearProblem)" begin
    @variables x(t) z(t)
    @mtkcompile sys = System(
        [D(x) ~ -x, 0 ~ atan(z - 3)], t; guesses = [z => 0.5]
    )
    iprob = InitializationProblem(sys, 0.0, [x => 1.0]; warn_initialize_determined = false)
    @test iprob isa SciMLBase.NonlinearProblem
    @test !(iprob isa SciMLBase.HomotopyProblem)
end

@testset "embedded ODEProblem stores a HomotopyProblem initializeprob" begin
    @variables x(t) y(t)
    @mtkcompile sys = System(
        [D(x) ~ -x, 0 ~ homotopy(atan(y - 3), y)], t; guesses = [y => 12.0]
    )
    odeprob = ODEProblem(sys, [x => 1.0], (0.0, 1.0))
    initprob = odeprob.f.initialization_data.initializeprob
    @test initprob isa SciMLBase.HomotopyProblem
    # the stored init problem solves to the true root by continuation
    isol = solve(initprob; abstol = 1.0e-12, reltol = 1.0e-12)
    @test SciMLBase.successful_retcode(isol)
    @test isol[y] ≈ 3.0 atol = 1.0e-6
end

@testset "homotopy init honors the requested in-place-ness" begin
    @variables x(t) y(t)
    @mtkcompile sys = System(
        [D(x) ~ -x, 0 ~ homotopy(atan(y - 3), y)], t; guesses = [y => 12.0]
    )
    # The init `op` is a varmap (never a StaticArray), so a bare `HomotopyProblem(isys, op)`
    # would always read as in-place; the branch must instead honor the caller's `{iip}`,
    # matching the sibling `NonlinearProblem{iip}` path.
    iip_prob = InitializationProblem{true}(sys, 0.0, [x => 1.0]; warn_initialize_determined = false)
    @test iip_prob isa SciMLBase.HomotopyProblem
    @test SciMLBase.isinplace(iip_prob)
    oop_prob = InitializationProblem{false}(sys, 0.0, [x => 1.0]; warn_initialize_determined = false)
    @test oop_prob isa SciMLBase.HomotopyProblem
    @test !SciMLBase.isinplace(oop_prob)
    # the out-of-place problem still solves by continuation
    sol = solve(oop_prob; abstol = 1.0e-12, reltol = 1.0e-12)
    @test SciMLBase.successful_retcode(sol)
    @test sol[y] ≈ 3.0 atol = 1.0e-6
end

@testset "get_initialization_problem_type selects by homotopy presence" begin
    @variables x(t) y(t)
    @mtkcompile hsys = System(
        [D(x) ~ -x, 0 ~ homotopy(atan(y - 3), y)], t; guesses = [y => 12.0]
    )
    # build the (square) init system the way InitializationProblem does and check the type
    isys = ModelingToolkitBase.generate_initializesystem(
        hsys; op = Dict(x => 1.0), guesses = Dict(y => 12.0)
    )
    isys = ModelingToolkitBase.mtkcompile(isys)
    @test has_any_homotopy(isys)
    @test get_initialization_problem_type(hsys, isys; warn_initialize_determined = false) ===
        SciMLBase.HomotopyProblem
    @test get_nonlinear_problem_type(isys) === SciMLBase.HomotopyProblem
end
