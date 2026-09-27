using ModelingToolkit, NonlinearSolve, SciMLBase, Test
using ModelingToolkitBase: t_nounits as t, SystemCompatibilityError

@testset "SCC two-blocks" begin
    @variables y(t) z(t)
    @parameters a b
    timedep = mtkcompile(
        System([0 ~ y^2 - a, 0 ~ z^2 - b], t, [y, z], [a, b]; name = :algebraic)
    )
    op = Dict(a => 9.0, b => 16.0)
    guesses = Dict(y => 2.5, z => 3.5)

    @testset "Multi-block time-dependent system fails check_compatibility" begin
        err = try
            SCCNonlinearProblem(timedep, op; combine_sccs = false, guesses)
            nothing
        catch e
            e
        end
        @test err isa SystemCompatibilityError
        @test occursin("mtkcompile(NonlinearSystem(sys))", sprint(showerror, err))
        @test !occursin("check_compatibility = false", sprint(showerror, err))
        @test_throws SystemCompatibilityError SCCNonlinearProblem(
            timedep, op; combine_sccs = true, guesses
        )
    end

    @testset "Single-block time-dependent falls back to NonlinearProblem" begin
        @variables w(t)
        @parameters c
        oneblock = mtkcompile(System([0 ~ w^2 - c], t, [w], [c]; name = :oneblock))
        @test length(ModelingToolkit.get_schedule(oneblock).var_sccs) == 1
        # Single-block path returns a `NonlinearProblem`, which takes unknowns in `op`
        # (same as constructing `NonlinearProblem` directly after `NonlinearSystem`).
        prob = SCCNonlinearProblem(oneblock, Dict(w => 2.5, c => 9.0))
        @test !(prob isa SciMLBase.SCCNonlinearProblem)
        sol = solve(prob; abstol = 1.0e-12, reltol = 1.0e-12)
        @test SciMLBase.successful_retcode(sol)
        @test sol[w] ≈ 3.0
    end

    @testset "After NonlinearSystem conversion, iip=$iip, combine=$combine" for iip in (
                true, false,
            ),
            combine in (true, false)
        sys = mtkcompile(NonlinearSystem(timedep))
        prob = SCCNonlinearProblem{iip}(
            sys, op; combine_sccs = combine, guesses
        )
        sol = solve(prob; abstol = 1.0e-12, reltol = 1.0e-12)
        @test SciMLBase.successful_retcode(sol)
        @test sol[y] ≈ 3.0
        @test sol[z] ≈ 4.0
    end
end
