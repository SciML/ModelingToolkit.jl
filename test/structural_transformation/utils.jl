using Test
using ModelingToolkit
using Graphs
using SparseArrays
using UnPack
using ModelingToolkit: t_nounits as t, D_nounits as D

# Define some variables
@parameters L g
@variables x(t) y(t) w(t) z(t) T(t)

# Simple pendulum in cartesian coordinates
eqs = [D(x) ~ w,
    D(y) ~ z,
    D(w) ~ T * x,
    D(z) ~ T * y - g,
    0 ~ x^2 + y^2 - L^2]
pendulum = ODESystem(eqs, t, [x, y, w, z, T], [L, g], name = :pendulum)
state = TearingState(pendulum)
StructuralTransformations.find_solvables!(state)
sss = state.structure
@unpack graph, solvable_graph, var_to_diff = sss
@test graph.fadjlist == [[1, 7], [2, 8], [3, 5, 9], [4, 6, 9], [5, 6]]
@test length(graph.badjlist) == 9
@test ne(graph) == nnz(incidence_matrix(graph)) == 12
@test nv(solvable_graph) == 9 + 5
let N = nothing
    @test var_to_diff == [N, N, N, N, 1, 2, 3, 4, N]
end

se = collect(StructuralTransformations.edges(graph))
@test se == mapreduce(vcat, enumerate(graph.fadjlist)) do (s, d)
    StructuralTransformations.BipartiteEdge.(s, d)
end

@testset "observed2graph handles unknowns inside callable parameters" begin
    @variables x(t) y(t)
    @parameters p(..)
    g, _ = ModelingToolkit.observed2graph([y ~ p(x), x ~ 0], [y, x])
    @test ModelingToolkit.𝑠neighbors(g, 1) == [2]
    @test ModelingToolkit.𝑑neighbors(g, 2) == [1]
end

@testset "array observed used unscalarized in another observed" begin
    @variables x(t) y(t)[1:2] z(t)[1:2]
    @parameters foo(::AbstractVector)[1:2]
    _tmp_fn(x) = 2x
    @mtkbuild sys = ODESystem(
        [D(x) ~ z[1] + z[2] + foo(z)[1], y[1] ~ 2t, y[2] ~ 3t, z ~ foo(y)], t)
    @test length(equations(sys)) == 1
    @test length(observed(sys)) == 7
    @test any(eq -> isequal(eq.lhs, y), observed(sys))
    @test any(eq -> isequal(eq.lhs, z), observed(sys))
    prob = ODEProblem(sys, [x => 1.0], (0.0, 1.0), [foo => _tmp_fn])
    @test_nowarn prob.f(prob.u0, prob.p, 0.0)

    isys = ModelingToolkit.generate_initializesystem(sys)
    @test length(unknowns(isys)) == 5
    @test length(equations(isys)) == 4
    @test !any(equations(isys)) do eq
        iscall(eq.rhs) && operation(eq.rhs) in [StructuralTransformations.getindex_wrapper,
            StructuralTransformations.change_origin]
    end
end

@testset "scalarized array observed calling same function multiple times" begin
    @variables x(t) y(t)[1:2]
    @parameters foo(::Real)[1:2]
    val = Ref(0)
    function _tmp_fn2(x)
        val[] += 1
        return [x, 2x]
    end
    @mtkbuild sys = ODESystem([D(x) ~ y[1] + y[2], y ~ foo(x)], t)
    @test length(equations(sys)) == 1
    @test length(observed(sys)) == 3
    prob = ODEProblem(sys, [x => 1.0], (0.0, 1.0), [foo => _tmp_fn2])
    @test_nowarn prob.f(prob.u0, prob.p, 0.0)
    @test val[] == 1

    isys = ModelingToolkit.generate_initializesystem(sys)
    @test length(unknowns(isys)) == 3
    @test length(equations(isys)) == 2
    @test !any(equations(isys)) do eq
        iscall(eq.rhs) && operation(eq.rhs) in [StructuralTransformations.getindex_wrapper,
            StructuralTransformations.change_origin]
    end

    @testset "CSE hack in equations(sys)" begin
        val[] = 0
        @variables z(t)[1:2]
        @mtkbuild sys = ODESystem(
            [D(y) ~ foo(x), D(x) ~ sum(y), zeros(2) ~ foo(prod(z))], t)
        @test length(equations(sys)) == 5
        @test length(observed(sys)) == 4
        prob = ODEProblem(
            sys, [y => ones(2), z => 2ones(2), x => 3.0], (0.0, 1.0), [foo => _tmp_fn2])
        @test_nowarn prob.f(prob.u0, prob.p, 0.0)
        @test val[] == 2

        isys = ModelingToolkit.generate_initializesystem(sys)
        @test length(unknowns(isys)) == 5
        @test length(equations(isys)) == 2
        @test !any(equations(isys)) do eq
            iscall(eq.rhs) &&
                operation(eq.rhs) in [StructuralTransformations.getindex_wrapper,
                    StructuralTransformations.change_origin]
        end
    end
end

@testset "array and cse hacks can be disabled" begin
    @testset "fully_determined = true" begin
        @variables x(t) y(t)[1:2] z(t)[1:2]
        @parameters foo(::AbstractVector)[1:2]
        _tmp_fn(x) = 2x
        @named sys = ODESystem(
            [D(x) ~ z[1] + z[2] + foo(z)[1], y[1] ~ 2t, y[2] ~ 3t, z ~ foo(y)], t)

        sys1 = structural_simplify(sys; cse_hack = false)
        @test length(observed(sys1)) == 6
        @test !any(observed(sys1)) do eq
            iscall(eq.rhs) &&
                operation(eq.rhs) == StructuralTransformations.getindex_wrapper
        end

        sys2 = structural_simplify(sys; array_hack = false)
        @test length(observed(sys2)) == 5
        @test !any(observed(sys2)) do eq
            iscall(eq.rhs) && operation(eq.rhs) == StructuralTransformations.change_origin
        end
    end

    @testset "fully_determined = false" begin
        @variables x(t) y(t)[1:2] z(t)[1:2] w(t)
        @parameters foo(::AbstractVector)[1:2]
        _tmp_fn(x) = 2x
        @named sys = ODESystem(
            [D(x) ~ z[1] + z[2] + foo(z)[1] + w, y[1] ~ 2t, y[2] ~ 3t, z ~ foo(y)], t)

        sys1 = structural_simplify(sys; cse_hack = false, fully_determined = false)
        @test length(observed(sys1)) == 6
        @test !any(observed(sys1)) do eq
            iscall(eq.rhs) &&
                operation(eq.rhs) == StructuralTransformations.getindex_wrapper
        end

        sys2 = structural_simplify(sys; array_hack = false, fully_determined = false)
        @test length(observed(sys2)) == 5
        @test !any(observed(sys2)) do eq
            iscall(eq.rhs) && operation(eq.rhs) == StructuralTransformations.change_origin
        end
    end
end

@testset "additional passes" begin
    @variables x(t) y(t)
    @named sys = ODESystem([D(x) ~ x, y ~ x + t], t)
    value = Ref(0)
    pass(sys; kwargs...) = (value[] += 1; return sys)
    structural_simplify(sys; additional_passes = [pass])
    @test value[] == 1
end
