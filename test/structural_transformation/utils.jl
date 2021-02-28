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
@test isequal(fullvars, [x, y, w, z, D(x), D(y), D(w), D(z), T])
@test graph.fadjlist == sort.([[5, 3], [6, 4], [7, 1, 9], [8, 2, 9], [2, 1]])
@test graph.badjlist == sort.([[3, 5], [4, 5], [1], [2], [1], [2], [3], [4], [3, 4]])
@test ne(graph) == nnz(incidence_matrix(graph)) == 12
@test nv(solvable_graph) == nv(solvable_graph) == 9 + 5
@test varassoc == [5, 6, 7, 8, 0, 0, 0, 0, 0]

se = collect(StructuralTransformations.𝑠edges(graph))
@test se == mapreduce(vcat, enumerate(graph.fadjlist)) do (s, d)
    StructuralTransformations.BipartiteEdge.(s, d)
end
de = collect(StructuralTransformations.𝑑edges(graph))
@test de == mapreduce(vcat, enumerate(graph.badjlist)) do (d, s)
    StructuralTransformations.BipartiteEdge.(s, d)
end
ae = collect(StructuralTransformations.edges(graph))
@test ae == vcat(se, de)
