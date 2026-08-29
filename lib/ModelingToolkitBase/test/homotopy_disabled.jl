using ModelingToolkitBase
using ModelingToolkitBase: homotopy, has_homotopy, has_any_homotopy, homotopy_enabled,
    strip_homotopy, get_nonlinear_problem_type, generate_rhs,
    t_nounits as t, D_nounits as D
using SciMLBase
using Symbolics
using Test

# `mtkcompile(sys; homotopy = false)` opts out of the continuation lowering: every
# `homotopy(actual, simplified)` node is replaced by `actual` before compilation, so the
# generated code only contains `actual`, problem construction never selects a
# `HomotopyProblem`, and systems derived from the compiled one inherit the setting.

@testset "strip_homotopy rewrites every field, recursively" begin
    @variables x(t) y(t) z(t)
    @named inner = System([0 ~ homotopy(y^2 - 2, y - 1.4)], t, [y], [])
    @named sys = System(
        [D(x) ~ homotopy(-x, -x + 1), z ~ homotopy(x^2, x)], t;
        systems = [inner], initialization_eqs = [homotopy(x^2, x) ~ 1],
        guesses = [x => homotopy(z + 1, z)]
    )
    ssys = strip_homotopy(sys)
    @test !any(eq -> has_homotopy(eq.lhs) || has_homotopy(eq.rhs), equations(ssys))
    @test !any(eq -> has_homotopy(eq.lhs) || has_homotopy(eq.rhs), initialization_equations(ssys))
    @test !any(has_homotopy, values(guesses(ssys)))
    @test isequal(equations(ssys)[1].rhs, -x)
    @test isequal(equations(ssys)[2].rhs, x^2)
    @test isequal(initialization_equations(ssys)[1].lhs, x^2)
    # the original is untouched
    @test any(eq -> has_homotopy(eq.rhs), equations(sys))
end

@testset "homotopy = false: generated code uses only the actual branch" begin
    @variables y
    @named sys = System([0 ~ homotopy(atan(y - 3), y)])
    csys = mtkcompile(sys; homotopy = false)
    @test !homotopy_enabled(csys)
    @test !has_any_homotopy(csys)

    expr = generate_rhs(csys; expression = Val{true})
    @test !occursin("homotopy", string(expr))
    @test occursin("atan", string(expr))

    @test get_nonlinear_problem_type(csys) === NonlinearProblem
    prob = SciMLBase.AbstractNonlinearProblem(csys, [y => 12.0])
    @test prob isa NonlinearProblem
    @test prob.f([12.0], prob.p) ≈ [atan(9.0)]
    @test_throws ArgumentError HomotopyProblem(csys, [y => 12.0])

    # default: nodes are kept and lowered
    dsys = mtkcompile(sys)
    @test homotopy_enabled(dsys)
    @test has_any_homotopy(dsys)
    @test occursin("homotopy", string(generate_rhs(dsys; expression = Val{true})))
    @test SciMLBase.AbstractNonlinearProblem(dsys, [y => 12.0]) isa HomotopyProblem
end

@testset "homotopy = false: initialization inherits the setting" begin
    @variables x(t) y(t)
    @named sys = System([D(x) ~ -x, 0 ~ homotopy(y^2 - x, y - 1)], t; guesses = [y => 1.5])
    csys = mtkcompile(sys; homotopy = false)
    prob = ODEProblem(csys, [x => 1.0], (0.0, 1.0))
    iprob = prob.f.initialization_data.initializeprob
    @test iprob isa NonlinearProblem
    @test !occursin("homotopy", string(generate_rhs(iprob.f.sys; expression = Val{true})))

    # initialization equations supplied at problem construction are stripped as well
    @named sys2 = System([D(x) ~ -x], t)
    csys2 = mtkcompile(sys2; homotopy = false)
    prob2 = ODEProblem(
        csys2, [], (0.0, 1.0);
        initialization_eqs = [homotopy(atan(x - 3) + x - 3, x - 3) ~ 0], guesses = [x => 1.0]
    )
    iprob2 = prob2.f.initialization_data.initializeprob
    @test iprob2 isa NonlinearProblem
    @test !has_any_homotopy(iprob2.f.sys)

    # with the default, the same problem routes to a HomotopyProblem
    prob3 = ODEProblem(
        mtkcompile(sys2), [], (0.0, 1.0);
        initialization_eqs = [homotopy(atan(x - 3) + x - 3, x - 3) ~ 0], guesses = [x => 1.0]
    )
    @test prob3.f.initialization_data.initializeprob isa HomotopyProblem
end
