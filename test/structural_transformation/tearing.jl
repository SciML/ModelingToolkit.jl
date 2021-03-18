using Test
using ModelingToolkit
using ModelingToolkit: Equation
using ModelingToolkit.StructuralTransformations: SystemStructure, find_solvables!
using NonlinearSolve
using LinearAlgebra
using UnPack

###
### Nonlinear system
###
@parameters t
@variables u1(t) u2(t) u3(t) u4(t) u5(t)
eqs = [
       0 ~ u1 - sin(u5),
       0 ~ u2 - cos(u1),
       0 ~ u3 - hypot(u1, u2),
       0 ~ u4 - hypot(u2, u3),
       0 ~ u5 - hypot(u4, u1),
]
sys = NonlinearSystem(eqs, [u1, u2, u3, u4, u5], [])
sys = initialize_system_structure(sys)
StructuralTransformations.find_solvables!(sys)
sss = structure(sys)
@unpack graph, solvable_graph, fullvars = sss

io = IOBuffer()
show(io, sss)
prt = String(take!(io))

if VERSION >= v"1.6"
@test prt == """xvars: Any[]
dxvars: Any[]
algvars: Any[u1(t), u2(t), u3(t), u4(t), u5(t)]
Incidence matrix:
 ×  ⋅  ⋅  ⋅  ×
 ×  ×  ⋅  ⋅  ⋅
 ×  ×  ×  ⋅  ⋅
 ⋅  ×  ×  ×  ⋅
 ×  ⋅  ⋅  ×  ×"""
end

# u1 = f1(u5)
# u2 = f2(u1)
# u3 = f3(u1, u2)
# u4 = f4(u2, u3)
# u5 = f5(u4, u1)
sys = initialize_system_structure(sys)
find_solvables!(sys)
sss = structure(sys)
@unpack graph, solvable_graph, assign, partitions = sss
@test graph.fadjlist == [[1, 5], [1, 2], [1, 2, 3], [2, 3, 4], [1, 4, 5]]
@test solvable_graph.fadjlist == map(x->[x], 1:5)

tornsys = tearing(sys)
sss = structure(tornsys)
@unpack graph, solvable_graph, assign, partitions = sss
@test graph.fadjlist == [[1]]
@test partitions == [([], [], [1], [1])]

# Before:
#      u1  u2  u3  u4  u5
# e1 [  1               1 ]
# e2 [  1   1             ]
# e3 [  1   1   1         ]
# e4 [      1   1   1     ]
# e5 [  1           1   1 ]
# solvable_graphs:
#      u1  u2  u3  u4  u5
# e1 [  1                 ]
# e2 [      1             ]
# e3 [          1         ]
# e4 [              1     ]
# e5 [                  1 ]
#
# Optimal:
#      u2  u3  u4  u5  | u1
# e2 [  1              |  1 ]
# e3 [  1   1          |  1 ]
# e4 [  1   1   1      |    ]
# e5 [          1   1  |  1 ]
# ---------------------|-----
# e1 [          1   1  |    ]
#
# Or:
#      u1  u2  u3  u4 | u5
# e1 [  1             |  1 ]
# e2 [  1   1         |    ]
# e3 [  1   1   1     |    ]
# e4 [      1   1   1 |    ]
# --------------------|-----
# e5 [  1           1 |  1 ]

sys = StructuralTransformations.tear_graph(StructuralTransformations.algebraic_equations_scc(sys))
sss = structure(sys)
@unpack partitions = sss
S = StructuralTransformations.reordered_matrix(sys, partitions)
@test S == [1 0 0 0 1
            1 1 0 0 0
            1 1 1 0 0
            0 1 1 1 0
            1 0 0 1 1]

# unknowns: u5
# u1 := sin(u5)
# u2 := cos(u1)
# u3 := hypot(u1, u2)
# u4 := hypot(u2, u3)
# solve for
# 0 = u5 - hypot(u1, u4)

# unknowns: u5
# solve for
# 0 = u5 - hypot(sin(u5), hypot(cos(sin(u5)), hypot(sin(u5), cos(sin(u5)))))
tornsys = tearing(sys)
@test isequal(equations(tornsys), [0 ~ u5 + (-1 * hypot(hypot(cos(sin(u5)), hypot(sin(u5), cos(sin(u5)))), sin(u5)))])
prob = NonlinearProblem(tornsys, ones(1))
sol = solve(prob, NewtonRaphson())
@test norm(prob.f(sol.u, sol.prob.p)) < 1e-10

###
### Simple test (edge case)
###
@parameters t
@variables x(t) y(t) z(t)
eqs = [
       0 ~ x - y,
       0 ~ z + y,
       0 ~ x + z,
      ]
nlsys = NonlinearSystem(eqs, [x, y, z], [])
newsys = tearing(nlsys)
@test equations(newsys) == [0 ~ z]
@test isequal(states(newsys), [z])

###
### DAE system
###
using ModelingToolkit, OrdinaryDiffEq, BenchmarkTools
@parameters t p
@variables x(t) y(t) z(t)
D = Differential(t)
eqs = [
       D(x) ~ z
       0 ~ x - y
       0 ~ sin(z) + y - p*t
      ]
daesys = ODESystem(eqs, t)
newdaesys = tearing(daesys)
@test equations(newdaesys) == [D(x) ~ z; 0 ~ x + sin(z) - p*t]
@test isequal(states(newdaesys), [x, z])
prob = ODAEProblem(newdaesys, [x=>1.0], (0, 1.0), [p=>0.2])
du = [0.0]; u = [1.0]; pr = 0.2; tt = 0.1
@test (@ballocated $(prob.f)($du, $u, $pr, $tt)) == 0
@test du ≈ [-asin(u[1] - pr * tt)] atol=1e-5

# test the initial guess is respected
infprob = ODAEProblem(tearing(ODESystem(eqs, t, defaults=Dict(z=>Inf))), [x=>1.0], (0, 1.0), [p=>0.2])
@test_throws DomainError infprob.f(du, u, pr, tt)

sol1 = solve(prob, Tsit5())
sol2 = solve(ODEProblem{false}(
                               (u,p,t) -> [-asin(u[1] - pr*t)],
                               [1.0],
                               (0, 1.0),
                               0.2,
                              ), Tsit5(), tstops=sol1.t, adaptive=false)
@test Array(sol1) ≈ Array(sol2) atol=1e-5

obs = build_observed_function(newdaesys, [z, y])
@test map(u -> u[2], obs.(sol1.u, pr, sol1.t)) == first.(sol1.u)
@test map(u -> sin(u[1]), obs.(sol1.u, pr, sol1.t)) + first.(sol1.u) ≈ pr[1]*sol1.t atol=1e-5

@test sol1[y, :] == sol1[x, :]
@test (@. sin(sol1[z, :]) + sol1[y, :]) ≈ pr * sol1.t atol=1e-5
