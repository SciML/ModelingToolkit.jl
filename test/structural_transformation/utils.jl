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
    @test length(observed(sys)) == 6
    @test any(eq -> isequal(eq.lhs, y), observed(sys))
    @test any(eq -> isequal(eq.lhs, z), observed(sys))
    prob = ODEProblem(sys, [x => 1.0], (0.0, 1.0), [foo => _tmp_fn])
    @test_nowarn prob.f(prob.u0, prob.p, 0.0)
end
