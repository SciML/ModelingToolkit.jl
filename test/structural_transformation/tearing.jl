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
@constants h = 1
@variables u1(t) u2(t) u3(t) u4(t) u5(t)
eqs = [
    0 ~ u1 - sin(u5) * h,
    0 ~ u2 - cos(u1),
    0 ~ u3 - hypot(u1, u2),
    0 ~ u4 - hypot(u2, u3),
    0 ~ u5 - hypot(u4, u1),
]
@named sys = NonlinearSystem(eqs, [u1, u2, u3, u4, u5], [])
state = TearingState(sys)
StructuralTransformations.find_solvables!(state)

io = IOBuffer()
show(io, MIME"text/plain"(), state.structure)
prt = String(take!(io))

if VERSION >= v"1.6"
    @test occursin("Incidence matrix:", prt)
    @test occursin("×", prt)
    @test occursin("⋅", prt)

    buff = IOBuffer()
    io = IOContext(buff, :mtk_limit => false)
    show(io, MIME"text/plain"(), state.structure)
    prt = String(take!(buff))
    @test occursin("SystemStructure", prt)
end

# u1 = f1(u5)
# u2 = f2(u1)
# u3 = f3(u1, u2)
# u4 = f4(u2, u3)
# u5 = f5(u4, u1)
state = TearingState(sys)
find_solvables!(state)
@unpack structure, fullvars = state
@unpack graph, solvable_graph = state.structure
int2var = Dict(eachindex(fullvars) .=> fullvars)
graph2vars(graph) = map(is -> Set(map(i -> int2var[i], is)), graph.fadjlist)
@test graph2vars(graph) == [Set([u1, u5])
       Set([u1, u2])
       Set([u1, u3, u2])
       Set([u4, u3, u2])
       Set([u4, u1, u5])]
@test graph2vars(solvable_graph) == [Set([u1])
                                     Set([u2])
                                     Set([u3])
                                     Set([u4])
                                     Set([u5])]

state = TearingState(tearing(sys))
let sss = state.structure
    @unpack graph = sss
    @test graph2vars(graph) == [Set([u1, u2, u5])]
end

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

let state = TearingState(sys)
    torn_matching = tearing(state)
    S = StructuralTransformations.reordered_matrix(sys, torn_matching)
    @test S == [1 0 0 0 1
                1 1 0 0 0
                1 1 1 0 0
                0 1 1 1 0
                1 0 0 1 1]
end

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
@test isequal(equations(tornsys), [0 ~ u5 - hypot(u4, u1)])
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
@named nlsys = NonlinearSystem(eqs, [x, y, z], [])
let (mm, _, _) = ModelingToolkit.aag_bareiss(nlsys)
    @test mm == [-1 1 0;
                 0 -1 -1;
                 0 0 0]
end

newsys = tearing(nlsys)
@test length(equations(newsys)) <= 1

###
### DAE system
###
using ModelingToolkit, OrdinaryDiffEq, BenchmarkTools
@parameters t p
@variables x(t) y(t) z(t)
D = Differential(t)
eqs = [D(x) ~ z * h
       0 ~ x - y
       0 ~ sin(z) + y - p * t]
@named daesys = ODESystem(eqs, t)
newdaesys = tearing(daesys)
@test equations(newdaesys) == [D(x) ~ h * z; 0 ~ y + sin(z) - p * t]
@test equations(tearing_substitution(newdaesys)) == [D(x) ~ h * z; 0 ~ x + sin(z) - p * t]
@test isequal(states(newdaesys), [x, z])
prob = ODAEProblem(newdaesys, [x => 1.0], (0, 1.0), [p => 0.2])
du = [0.0];
u = [1.0];
pr = 0.2;
tt = 0.1;
@test_skip (@ballocated $(prob.f)($du, $u, $pr, $tt)) == 0
prob.f(du, u, pr, tt)
@test du≈[-asin(u[1] - pr * tt)] atol=1e-5

# test the initial guess is respected
@named sys = ODESystem(eqs, t, defaults = Dict(z => Inf))
infprob = ODAEProblem(tearing(sys), [x => 1.0], (0, 1.0), [p => 0.2])
@test_throws DomainError infprob.f(du, u, pr, tt)

sol1 = solve(prob, Tsit5())
sol2 = solve(ODEProblem{false}((u, p, t) -> [-asin(u[1] - pr * t)],
                               [1.0],
                               (0, 1.0),
                               0.2), Tsit5(), tstops = sol1.t, adaptive = false)
@test Array(sol1)≈Array(sol2) atol=1e-5

@test sol1[x] == first.(sol1.u)
@test sol1[y] == first.(sol1.u)
@test sin.(sol1[z]) .+ sol1[y]≈pr[1] * sol1.t atol=1e-5
@test sol1[sin(z) + y]≈sin.(sol1[z]) .+ sol1[y] rtol=1e-12

@test sol1[y, :] == sol1[x, :]
@test (@. sin(sol1[z, :]) + sol1[y, :])≈pr * sol1.t atol=1e-5

# 1426
function Translational_Mass(; name, m = 1.0)
    sts = @variables s(t) v(t) a(t)
    ps = @parameters m = m
    D = Differential(t)
    eqs = [D(s) ~ v
           D(v) ~ a
           m * a ~ 0.0]
    ODESystem(eqs, t, sts, ps; name = name)
end

m = 1.0
@named mass = Translational_Mass(m = m)

ms_eqs = []

@named _ms_model = ODESystem(ms_eqs, t)
@named ms_model = compose(_ms_model,
                          [mass])

calculate_jacobian(ms_model)
calculate_tgrad(ms_model)

# Mass starts with velocity = 1
u0 = [mass.s => 0.0
      mass.v => 1.0]

sys = structural_simplify(ms_model)
@test ModelingToolkit.get_jac(sys)[] === ModelingToolkit.EMPTY_JAC
@test ModelingToolkit.get_tgrad(sys)[] === ModelingToolkit.EMPTY_TGRAD
prob_complex = ODAEProblem(sys, u0, (0, 1.0))
sol = solve(prob_complex, Tsit5())
@test all(sol[mass.v] .== 1)
