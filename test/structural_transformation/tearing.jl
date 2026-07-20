using Test
using ModelingToolkit
using ModelingToolkit: Equation, observed
using ModelingToolkit.StructuralTransformations: SystemStructure
using NonlinearSolve
using LinearAlgebra
using UnPack
using SymbolicIndexingInterface
using ModelingToolkit: t_nounits as t, D_nounits as D
import StateSelection
import SymbolicUtils as SU
using ForwardDiff

###
### Nonlinear system
###
@constants h = 1
@variables u1(t) u2(t) u3(t) u4(t) u5(t) [state_priority = 10]
eqs = [
    0 ~ u1 - sin(u5) * h,
    0 ~ u2 - cos(u1),
    0 ~ u3 - hypot(u1, u2),
    0 ~ u4 - hypot(u2, u3),
    0 ~ u5 - hypot(u4, u1),
]
@named sys = System(eqs, [u1, u2, u3, u4, u5], [h])
state = TearingState(sys)
StateSelection.find_solvables!(state)

io = IOBuffer()
show(io, MIME"text/plain"(), state.structure)
prt = String(take!(io))

@test occursin("Incidence matrix:", prt)
@test occursin("×", prt)
@test occursin("⋅", prt)

buff = IOBuffer()
io = IOContext(buff, :mtk_limit => false)
show(io, MIME"text/plain"(), state.structure)
prt = String(take!(buff))
@test occursin("SystemStructure", prt)

# u1 = f1(u5)
# u2 = f2(u1)
# u3 = f3(u1, u2)
# u4 = f4(u2, u3)
# u5 = f5(u4, u1)
state = TearingState(sys)
StateSelection.find_solvables!(state)
@unpack structure, fullvars = state
@unpack graph, solvable_graph = state.structure
int2var = Dict(eachindex(fullvars) .=> fullvars)
graph2vars(graph) = map(is -> Set(map(i -> int2var[i], is)), graph.fadjlist)
# @test graph2vars(graph) == [
#         Set([u4, u3, u2])
#         Set([u1, u5])
#         Set([u1, u2])
#         Set([u1, u2, u3])
#         Set([u4, u1, u5])
#     ]
# @test graph2vars(solvable_graph) == [Set([u4])
#                                      Set([u1])
#                                      Set([u2])
#                                      Set([u3])
#                                      Set([u5])]

newsys = tearing(sys)
@test length(equations(newsys)) == 1
vars = SU.search_variables(equations(newsys))
@test issetequal(vars, [u1, u4, u5]) || issetequal(vars, [h, u1, u5])

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

# Variable and equation ordering changes with version
if v"1.12" <= VERSION < v"1.13"
    let state = TearingState(sys)
        result, = tearing(state)
        S = StructuralTransformations.reordered_matrix(sys, result.var_eq_matching)
        # Tearing `u5` or tearing `u1` (the two arrangements diagrammed above) are both
        # valid maximal matchings; which one the heuristic picks depends on the
        # lexicographic equation sort in `TearingState`, whose strings depend on the
        # symbolic backend's term ordering (e.g. it flipped with a SymbolicUtils hashing
        # change, see #4619). Accept either, like the `vars` test above does.
        S_tear_u5 = [
            1 0 0 0 1
            1 1 0 0 0
            1 1 1 0 0
            0 1 1 1 0
            1 0 0 1 1
        ]
        S_tear_u1 = [
            1 0 0 1 1
            0 1 0 0 1
            0 1 1 0 1
            0 1 1 1 0
            1 0 0 0 1
        ]
        @test S == S_tear_u5 || S == S_tear_u1
    end
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
tornsys = complete(tearing(sys))
@test isequal(equations(tornsys), [0 ~ u5 - hypot(u4, u1)]) || isequal(equations(tornsys), [0 ~ u1 - h * sin(u5)])
prob = NonlinearProblem(tornsys, [u5 => 1.0, u1 => 1.0])
sol = solve(prob, NewtonRaphson())
@test norm(prob.f(sol.u, sol.prob.p)) < 1.0e-10

###
### Simple test (edge case)
###
@variables x(t) y(t) z(t)
eqs = [
    0 ~ x - y,
    0 ~ z + y,
    0 ~ x + z,
]
@named nlsys = System(eqs, [x, y, z], [])

newsys = tearing(nlsys)
@test length(equations(newsys)) <= 1

###
### DAE system
###
using ModelingToolkit, OrdinaryDiffEq, BenchmarkTools
@parameters p
@variables x(t) [state_priority = 1] y(t) z(t)
eqs = [
    D(x) ~ z * h
    0 ~ x - y
    0 ~ sin(z) + y - p * t
]
@named daesys = System(eqs, t)
newdaesys = mtkcompile(daesys)
@test issetequal(equations(newdaesys), [D(x) ~ h * z; 0 ~ x + sin(z) - p * t])
@test issetequal(
    equations(tearing_substitution(newdaesys)), [D(x) ~ h * z; 0 ~ x + sin(z) - p * t]
)
@test issetequal(unknowns(newdaesys), [x, z])
prob = ODEProblem(newdaesys, [x => 1.0, z => -0.5π, p => 0.2], (0, 1.0))
du = [0.0, 0.0];
u = [1.0, -0.5π];
pr = prob.p;
tt = 0.1;
@test (@ballocated $(prob.f)($du, $u, $pr, $tt)) == 0
prob.f(du, u, pr, tt)
xgetter = getsym(prob, x)
zgetter = getsym(prob, z)
@test xgetter(du) ≈ zgetter(u) atol = 1.0e-5
@test zgetter(du) ≈ xgetter(u) + sin(zgetter(u)) - prob.ps[p] * tt atol = 1.0e-5

# test the initial guess is respected
@named sys = System(eqs, t, initial_conditions = Dict(z => NaN))
infprob = ODEProblem(mtkcompile(sys), [x => 1.0, p => 0.2], (0, 1.0))
infprob.f(du, infprob.u0, pr, tt)
@test any(isnan, du)

# 1426
function Translational_Mass(; name, m = 1.0)
    sts = @variables s(t) v(t) a(t)
    ps = @parameters m = m
    D = Differential(t)
    eqs = [
        D(s) ~ v
        D(v) ~ a
        m * a ~ 0.0
    ]
    return System(eqs, t, sts, ps; name = name)
end

m = 1.0
@named mass = Translational_Mass(m = m)

ms_eqs = Equation[]

@named _ms_model = System(ms_eqs, t)
@named ms_model = compose(
    _ms_model,
    [mass]
)

calculate_jacobian(ms_model)
calculate_tgrad(ms_model)

# Mass starts with velocity = 1
u0 = [
    mass.s => 0.0
    mass.v => 1.0
]

sys = mtkcompile(ms_model)
# @test ModelingToolkit.get_jac(sys)[] === ModelingToolkit.EMPTY_JAC
# @test ModelingToolkit.get_tgrad(sys)[] === ModelingToolkit.EMPTY_TGRAD
prob_complex = ODEProblem(sys, u0, (0, 1.0))
sol = solve(prob_complex, Tsit5())
@test all(sol[mass.v] .== 1)

@testset "Inline linear SCCs" begin
    reassemble_alg0 = StructuralTransformations.DefaultReassembleAlgorithm(; inline_linear_sccs = false)
    reassemble_alg1 = StructuralTransformations.DefaultReassembleAlgorithm(; inline_linear_sccs = true, analytical_linear_scc_limit = 1)
    reassemble_alg2 = StructuralTransformations.DefaultReassembleAlgorithm(; inline_linear_sccs = true, analytical_linear_scc_limit = 5)
    @variables x(t)[1:3] y(t)[1:3]
    @parameters p[1:3, 1:3]
    @mtkcompile sys1 = System([D(x) ~ x, p * y ~ x], t) reassemble_alg = reassemble_alg0
    @mtkcompile sys2 = System([D(x) ~ x, p * y ~ x], t) reassemble_alg = reassemble_alg1
    @mtkcompile sys3 = System([D(x) ~ x, p * y ~ x], t) reassemble_alg = reassemble_alg2

    @test length(equations(sys1)) > 3
    @test length(equations(sys2)) == 3
    @test length(equations(sys3)) == 3

    idx = findfirst(observed(sys2)) do eq
        arr, isarr = ModelingToolkit.MTKBase.split_indexed_var(eq.rhs)
        isarr && iscall(arr) && operation(arr) === (\)
    end
    @test idx !== nothing
    # This one is analytically solved
    idx = findfirst(observed(sys3)) do eq
        arr, isarr = ModelingToolkit.MTKBase.split_indexed_var(eq.rhs)
        isarr && iscall(arr) && operation(arr) === (\)
    end
    @test idx === nothing

    pval = rand(3, 3)
    prob1 = ODEProblem(sys1, [x => ones(3), p => pval], (0.0, 10.0); guesses = [y => ones(3)])
    prob2 = ODEProblem(sys2, [x => ones(3), p => pval], (0.0, 10.0))
    prob3 = ODEProblem(sys3, [x => ones(3), p => pval], (0.0, 10.0))

    sol1 = solve(prob1, Rodas5P(); abstol = 1.0e-8, reltol = 1.0e-8)
    sol2 = solve(prob2, Tsit5(), abstol = 1.0e-8, reltol = 1.0e-8)
    sol3 = solve(prob3, Tsit5(), abstol = 1.0e-8, reltol = 1.0e-8)

    @test SciMLBase.successful_retcode(sol1)
    @test SciMLBase.successful_retcode(sol2)
    @test SciMLBase.successful_retcode(sol3)

    @test sol2(sol1.t; idxs = unknowns(sys1)).u ≈ sol1.u rtol = 1.0e-6
    @test sol3(sol1.t; idxs = unknowns(sys1)).u ≈ sol1.u rtol = 1.0e-6
end

@testset "`Initial` parameters are added for observed variables solved by inline linear SCCS" begin
    reassemble_alg = StructuralTransformations.DefaultReassembleAlgorithm(; inline_linear_sccs = true, analytical_linear_scc_limit = 1)
    @variables x(t)[1:3] y(t)[1:3]
    @parameters p[1:3, 1:3]
    @mtkcompile sys = System([D(x) ~ x, p * y ~ x], t) reassemble_alg = reassemble_alg
    @test Initial(y) in Set(ModelingToolkit.get_ps(sys))
end

@testset "AD through inline linear SCCs works" begin
    reassemble_alg = StructuralTransformations.DefaultReassembleAlgorithm(; inline_linear_sccs = true, analytical_linear_scc_limit = 1)
    @variables x(t)[1:3] y(t)[1:3]
    @parameters p[1:3, 1:3]
    @mtkcompile sys = System([D(x) ~ x, p * y ~ x], t) reassemble_alg = reassemble_alg
    prob = ODEProblem(sys, [x => [1.0, 2.3, 5.7], p => rand(3, 3)], (0.0, 10.0))
    @assert prob.p.nonnumeric[1] isa
        Vector{ModelingToolkitBase.DiffCacheAllocatorAPIWrapper{Float64}}
    @assert SciMLBase.has_initializeprob(prob.f)

    setter = setsym_oop(prob, [p[1, 1]])

    function loss(x)
        new_u0, new_p = setter(prob, x)
        new_prob = remake(prob; p = new_p)
        sol = solve(new_prob, Tsit5(); abstol = 1.0e-8, reltol = 1.0e-8)
        return sol[sys.y[1]][end]
    end

    # Primal works: returns ~0.993
    @test_nowarn loss([1.0])

    @test_nowarn ForwardDiff.gradient(loss, [1.0])
end
