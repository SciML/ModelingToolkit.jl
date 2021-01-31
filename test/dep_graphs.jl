using Test
using ModelingToolkit, LightGraphs

import ModelingToolkit: value

# use a ReactionSystem to generate systems for testing
@parameters k1 k2 t
@variables S(t) I(t) R(t)

rxs = [Reaction(k1, nothing, [S]),
       Reaction(k1, [S], nothing),
       Reaction(k2, [S,I], [I], [1,1], [2]),
       Reaction(k2, [S,R], [S], [2,1], [2]),
       Reaction(k1*I, nothing, [R]),
       Reaction(k1*k2/(1+t), [S], [R])]
rs = ReactionSystem(rxs, t, [S,I,R], [k1,k2])


#################################
#  testing for Jumps / all dgs
#################################
js = convert(JumpSystem, rs)
S  = value(S); I = value(I); R = value(R)
k1 = value(k1); k2 = value(k2)
# eq to vars they depend on
eq_sdeps  = [Variable[], [S], [S,I], [S,R], [I], [S]]
eq_sidepsf = [Int[], [1], [1,2], [1,3], [2], [1]]
eq_sidepsb = [[2,3,4,6], [3,5],[4]]
deps = equation_dependencies(js)
@test all(i -> isequal(Set(eq_sdeps[i]),Set(deps[i])), 1:length(rxs))
depsbg = asgraph(js)
@test depsbg.fadjlist == eq_sidepsf
@test depsbg.badjlist == eq_sidepsb

# eq to params they depend on
eq_pdeps    = [[k1],[k1],[k2],[k2],[k1],[k1,k2]]
eq_pidepsf = [[1],[1],[2],[2],[1],[1,2]]
eq_pidepsb = [[1,2,5,6],[3,4,6]]
deps = equation_dependencies(js, variables=parameters(js))
@test all(i -> isequal(Set(eq_pdeps[i]),Set(deps[i])), 1:length(rxs))
depsbg2 = asgraph(js, variables=parameters(js))
@test depsbg2.fadjlist == eq_pidepsf
@test depsbg2.badjlist == eq_pidepsb

# var to eqs that modify them
s_eqdepsf = [[1,2,3,6],[3],[4,5,6]]
s_eqdepsb = [[1],[1],[1,2],[3],[3],[1,3]]
ne        = 8
bg        = BipartiteGraph(ne, s_eqdepsf, s_eqdepsb)
deps2     = variable_dependencies(js)
@test isequal(bg,deps2)

# eq to eqs that depend on them
eq_eqdeps = [[2,3,4,6],[2,3,4,6],[2,3,4,5,6],[4],[4],[2,3,4,6]]
dg = SimpleDiGraph(6)
for (eqidx,eqdeps) in enumerate(eq_eqdeps)
    for eqdepidx in eqdeps
       add_edge!(dg, eqidx, eqdepidx)
    end
end
dg3 = eqeq_dependencies(depsbg,deps2)
@test dg == dg3

# var to vars that depend on them
var_vardeps = [[1,2,3],[1,2,3],[3]]
ne = 7
dg = SimpleDiGraph(3)
for (vidx,vdeps) in enumerate(var_vardeps)
    for vdepidx in vdeps
       add_edge!(dg, vidx, vdepidx)
    end
end
dg4 = varvar_dependencies(depsbg,deps2)
@test dg == dg4

#####################################
#       testing for ODE/SDEs
#####################################
os = convert(ODESystem, rs)
deps = equation_dependencies(os)
eq_sdeps  = [[S,I], [S,I], [S,I,R]]
@test all(i -> isequal(Set(eq_sdeps[i]),Set(deps[i])), 1:length(deps))

sdes = convert(SDESystem, rs)
deps = equation_dependencies(sdes)
@test all(i -> isequal(Set(eq_sdeps[i]),Set(deps[i])), 1:length(deps))

deps = variable_dependencies(os)
s_eqdeps = [[1],[2],[3]]
@test deps.fadjlist == s_eqdeps

#####################################
#       testing for nonlin sys
#####################################
@variables x y z
@parameters σ ρ β

eqs = [0 ~ σ*(y-x),
       0 ~ ρ-y,
       0 ~ y - β*z]
ns = NonlinearSystem(eqs, [x,y,z],[σ,ρ,β])
deps = equation_dependencies(ns)
eq_sdeps = [[x,y],[y],[y,z]]
@test all(i -> isequal(Set(deps[i]),Set(value.(eq_sdeps[i]))), 1:length(deps))

using SparseArrays
using ModelingToolkit
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
sss = structure(sys)
@unpack graph, fullvars, varassoc = sss
@test isequal(fullvars, [x, y, w, z, D(x), D(y), D(w), D(z), T])
@test graph.fadjlist == sort.([[5, 3], [6, 4], [7, 1, 9], [8, 2, 9], [2, 1]])
@test graph.badjlist == sort.([[3, 5], [4, 5], [1], [2], [1], [2], [3], [4], [3, 4]])
@test LightGraphs.ne(graph) == nnz(incidence_matrix(graph)) == 12
@test varassoc == [5, 6, 7, 8, 0, 0, 0, 0, 0]

se = collect(ModelingToolkit.𝑠edges(graph))
@test se == mapreduce(vcat, enumerate(graph.fadjlist)) do (s, d)
    ModelingToolkit.BipartiteEdge.(s, d)
end
de = collect(ModelingToolkit.𝑑edges(graph))
@test de == mapreduce(vcat, enumerate(graph.badjlist)) do (d, s)
    ModelingToolkit.BipartiteEdge.(s, d)
end
ae = collect(ModelingToolkit.edges(graph))
@test ae == vcat(se, de)
