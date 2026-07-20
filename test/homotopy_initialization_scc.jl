using ModelingToolkit
using ModelingToolkit: homotopy, t_nounits as t, D_nounits as D
using SciMLBase
using NonlinearSolve   # continuation solvers for the HomotopyProblem
using Test

# The full-MTK `get_initialization_problem_type` selects `SCCNonlinearProblem` for a
# fully determined, split initialization system that tears into more than one block.
# Modelica's `homotopy(actual, simplified)` shares one continuation parameter λ across the
# whole model (spec 3.7.4.2), so a `homotopy`-carrying init system must instead be solved
# as a single `HomotopyProblem` rather than torn into SCC blocks (each of which would
# sweep its own λ with per-block Newton solves).
#
# The fixture is a chain `a -> b -> c` of algebraic equations, which tears into three
# sequential SCC blocks — a single-equation init system degrades to a plain
# `NonlinearProblem` inside the `SCCNonlinearProblem` constructor and would not exercise
# this path at all.

@testset "homotopy init system bypasses SCC and builds a HomotopyProblem" begin
    @variables x(t) a(t) b(t) c(t)
    @mtkcompile sys = System(
        [D(x) ~ -x, 0 ~ homotopy(a^2 - 2, a - 1.414), 0 ~ b^2 - a, 0 ~ c^2 - b],
        t; guesses = [a => 1.5, b => 1.2, c => 1.1]
    )
    @test ModelingToolkit.is_split(sys)

    iprob = ModelingToolkit.InitializationProblem(
        sys, 0.0, [x => 1.0]; warn_initialize_determined = false
    )
    @test iprob isa SciMLBase.HomotopyProblem
    @test !(iprob isa SciMLBase.SCCNonlinearProblem)

    # the whole (3-block) init system is solved by continuation, to the exact roots
    sol = solve(iprob; abstol = 1.0e-12, reltol = 1.0e-12)
    @test SciMLBase.successful_retcode(sol)
    @test sol[a] ≈ sqrt(2) atol = 1.0e-8
    @test sol[b] ≈ sqrt(sqrt(2)) atol = 1.0e-8
    @test sol[c] ≈ sqrt(sqrt(sqrt(2))) atol = 1.0e-8
end

@testset "non-homotopy init system still uses SCCNonlinearProblem" begin
    # regression guard: the SCC path must be untouched for systems without `homotopy`
    @variables x(t) a(t) b(t) c(t)
    @mtkcompile sys = System(
        [D(x) ~ -x, 0 ~ a^2 - 2, 0 ~ b^2 - a, 0 ~ c^2 - b],
        t; guesses = [a => 1.5, b => 1.2, c => 1.1]
    )
    iprob = ModelingToolkit.InitializationProblem(
        sys, 0.0, [x => 1.0]; warn_initialize_determined = false
    )
    @test iprob isa SciMLBase.SCCNonlinearProblem
    @test !(iprob isa SciMLBase.HomotopyProblem)
end

@testset "use_scc = false still routes homotopy to HomotopyProblem" begin
    @variables x(t) a(t) b(t) c(t)
    @mtkcompile sys = System(
        [D(x) ~ -x, 0 ~ homotopy(a^2 - 2, a - 1.414), 0 ~ b^2 - a, 0 ~ c^2 - b],
        t; guesses = [a => 1.5, b => 1.2, c => 1.1]
    )
    iprob = ModelingToolkit.InitializationProblem(
        sys, 0.0, [x => 1.0]; warn_initialize_determined = false, use_scc = false
    )
    @test iprob isa SciMLBase.HomotopyProblem
end
