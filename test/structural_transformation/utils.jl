using Test
using ModelingToolkit
using Graphs
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
@test isequal(fullvars, [D(x), D(y), D(w), D(z), x, y, w, z, T])
@test graph.fadjlist == [[1, 7], [2, 8], [3, 5, 9], [4, 6, 9], [5, 6]]
@test graph.badjlist == 9 == length(fullvars)
@test ne(graph) == nnz(incidence_matrix(graph)) == 12
@test nv(solvable_graph) == 9 + 5
@test varassoc == [0, 0, 0, 0, 1, 2, 3, 4, 0]

se = collect(StructuralTransformations.edges(graph))
@test se == mapreduce(vcat, enumerate(graph.fadjlist)) do (s, d)
    StructuralTransformations.BipartiteEdge.(s, d)
end
