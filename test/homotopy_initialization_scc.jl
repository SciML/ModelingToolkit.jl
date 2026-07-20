using ModelingToolkit
using ModelingToolkit: homotopy, t_nounits as t, D_nounits as D
using SciMLBase
using NonlinearSolve   # continuation solvers for the HomotopyProblem blocks
using Test

# Modelica `homotopy(actual, simplified)` nodes are handled PER SCC: the initialization
# system still tears into `SCCNonlinearProblem` blocks, but a block whose equations carry
# `homotopy` nodes is built as a `SciMLBase.HomotopyProblem` (λ-swept residual) and solved
# by continuation, while the remaining blocks keep their plain Newton solves. Blocks are
# solved in dependency order, so each reaches `λ = 1` before the next begins.
#
# The fixture is a chain `a -> b -> c` of algebraic equations, which tears into three
# sequential SCC blocks — a single-equation init system degrades to a plain
# `NonlinearProblem` inside the `SCCNonlinearProblem` constructor and would not exercise
# this path at all. Only the `a` block carries a `homotopy` node.

@testset "homotopy is applied per SCC block" begin
    @variables x(t) a(t) b(t) c(t)
    @mtkcompile sys = System(
        [D(x) ~ -x, 0 ~ homotopy(a^2 - 2, a - 1.414), 0 ~ b^2 - a, 0 ~ c^2 - b],
        t; guesses = [a => 1.5, b => 1.2, c => 1.1]
    )
    @test ModelingToolkit.is_split(sys)

    iprob = ModelingToolkit.InitializationProblem(
        sys, 0.0, [x => 1.0]; warn_initialize_determined = false
    )
    # the init system is still torn into SCC blocks ...
    @test iprob isa SciMLBase.SCCNonlinearProblem
    # ... and exactly the `homotopy`-carrying block is a HomotopyProblem
    nhomotopy = count(p -> p isa SciMLBase.HomotopyProblem, iprob.probs)
    @test nhomotopy == 1
    @test count(p -> p isa SciMLBase.NonlinearProblem, iprob.probs) ==
        length(iprob.probs) - 1

    # the λ-swept block residual genuinely depends on λ
    hblock = iprob.probs[findfirst(p -> p isa SciMLBase.HomotopyProblem, iprob.probs)]
    u = copy(SciMLBase.state_values(hblock))
    pp = SciMLBase.parameter_values(hblock)
    if SciMLBase.isinplace(hblock)
        r0 = similar(u)
        r1 = similar(u)
        hblock.f(r0, u, pp, 0.0)
        hblock.f(r1, u, pp, 1.0)
        @test !(r0 ≈ r1)
    else
        @test !(hblock.f(u, pp, 0.0) ≈ hblock.f(u, pp, 1.0))
    end

    # solving the SCC problem walks the homotopy block by continuation and the rest by
    # Newton, landing on the exact roots
    sol = solve(iprob; abstol = 1.0e-12, reltol = 1.0e-12)
    @test SciMLBase.successful_retcode(sol)
    @test sol[a] ≈ sqrt(2) atol = 1.0e-8
    @test sol[b] ≈ sqrt(sqrt(2)) atol = 1.0e-8
    @test sol[c] ≈ sqrt(sqrt(sqrt(2))) atol = 1.0e-8
end

@testset "non-homotopy init system has no HomotopyProblem blocks" begin
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
    @test !any(p -> p isa SciMLBase.HomotopyProblem, iprob.probs)
end

@testset "use_scc = false falls back to a whole-system HomotopyProblem" begin
    # with no SCC decomposition to hang per-block continuation off, the whole init system
    # is solved as a single HomotopyProblem
    @variables x(t) a(t) b(t) c(t)
    @mtkcompile sys = System(
        [D(x) ~ -x, 0 ~ homotopy(a^2 - 2, a - 1.414), 0 ~ b^2 - a, 0 ~ c^2 - b],
        t; guesses = [a => 1.5, b => 1.2, c => 1.1]
    )
    iprob = ModelingToolkit.InitializationProblem(
        sys, 0.0, [x => 1.0]; warn_initialize_determined = false, use_scc = false
    )
    @test iprob isa SciMLBase.HomotopyProblem
    sol = solve(iprob; abstol = 1.0e-12, reltol = 1.0e-12)
    @test SciMLBase.successful_retcode(sol)
    @test sol[a] ≈ sqrt(2) atol = 1.0e-8
end
