using Test
using ModelingToolkit
using LightGraphs
using SparseArrays
using UnPack

# Define some variables
@parameters t L g
@variables x(t) y(t) w(t) z(t) T(t)
D = Differential(t)

# Simple pendulum in cartesian coordinates
eqs = [D(x) ~ w,
       D(y) ~ z,
       D(w) ~ T*x,
       D(z) ~ T*y - g,
       0 ~ x^2 + y^2 - L^2]
pendulum = ODESystem(eqs, t, [x, y, w, z, T], [L, g], name=:pendulum)
sys = initialize_system_structure(pendulum)
StructuralTransformations.find_solvables!(sys)
sss = structure(sys)
@unpack graph, solvable_graph, fullvars, varassoc = sss
@test isequal(fullvars, [w, D(x), x, z, D(y), y, T, D(w), D(z)])
@test graph.fadjlist == sort.([[1, 2], [4, 5], [3, 7, 8], [6, 7, 9], [3, 6]])
@test graph.badjlist == 9 == length(fullvars)
@test ne(graph) == nnz(incidence_matrix(graph)) == 12
@test nv(solvable_graph) == 9 + 5
@test varassoc == [8, 0, 2, 9, 0, 5, 0, 0, 0]

se = collect(StructuralTransformations.𝑠edges(graph))
@test se == mapreduce(vcat, enumerate(graph.fadjlist)) do (s, d)
    StructuralTransformations.BipartiteEdge.(s, d)
end
@test_throws ArgumentError collect(StructuralTransformations.𝑑edges(graph))
@test_throws ArgumentError collect(StructuralTransformations.edges(graph))
