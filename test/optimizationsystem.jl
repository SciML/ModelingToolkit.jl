using ModelingToolkit, SparseArrays, Test, Optimization, OptimizationOptimJL,
      OptimizationMOI, Ipopt, AmplNLWriter, Ipopt_jll
using ModelingToolkit: get_metadata

@variables x y
@parameters a b
loss = (a - x)^2 + b * (y - x^2)^2
sys1 = OptimizationSystem(loss, [x, y], [a, b], name = :sys1)

cons2 = [x^2 + y^2 ~ 0, y * sin(x) - x ~ 0]
sys2 = OptimizationSystem(loss, [x, y], [a, b], name = :sys2, constraints = cons2)

@variables z
@parameters β
loss2 = sys1.x - sys2.y + z * β
combinedsys = OptimizationSystem(loss2, [z], [β], systems = [sys1, sys2],
                                 name = :combinedsys)

equations(combinedsys)
states(combinedsys)
parameters(combinedsys)

calculate_gradient(combinedsys)
calculate_hessian(combinedsys)
generate_function(combinedsys)
generate_gradient(combinedsys)
generate_hessian(combinedsys)
hess_sparsity = ModelingToolkit.hessian_sparsity(sys1)
sparse_prob = OptimizationProblem(sys1, [x, y], [a, b], grad = true, sparse = true)
@test sparse_prob.f.hess_prototype.rowval == hess_sparsity.rowval
@test sparse_prob.f.hess_prototype.colptr == hess_sparsity.colptr

u0 = [sys1.x => 1.0
      sys1.y => 2.0
      sys2.x => 3.0
      sys2.y => 4.0
      z => 5.0]
p = [sys1.a => 6.0
     sys1.b => 7.0
     sys2.a => 8.0
     sys2.b => 9.0
     β => 10.0]

prob = OptimizationProblem(combinedsys, u0, p, grad = true)
@test prob.f.sys === combinedsys
sol = solve(prob, NelderMead())
@test sol.minimum < -1e5

prob2 = remake(prob, u0 = sol.minimizer)
sol = solve(prob, BFGS(initial_stepnorm = 0.0001), allow_f_increases = true)
@test sol.minimum < -1e8
sol = solve(prob2, BFGS(initial_stepnorm = 0.0001), allow_f_increases = true)
@test sol.minimum < -1e9

#inequality constraint, the bounds for constraints lcons !== ucons
prob = OptimizationProblem(sys2, [x => 0.0, y => 0.0], [a => 1.0, b => 100.0],
                           lcons = [-1.0, -1.0], ucons = [500.0, 500.0], grad = true,
                           hess = true)
@test prob.f.sys === sys2
sol = solve(prob, IPNewton(), allow_f_increases = true)
@test sol.minimum < 1.0
sol = solve(prob, Ipopt.Optimizer())
@test sol.minimum < 1.0
sol = solve(prob, AmplNLWriter.Optimizer(Ipopt_jll.amplexe))
@test sol.minimum < 1.0

#equality constraint, lcons == ucons
cons2 = [0.0 ~ x^2 + y^2]
out = zeros(1)
sys2 = OptimizationSystem(loss, [x, y], [a, b], name = :sys2, constraints = cons2)
prob = OptimizationProblem(sys2, [x => 0.0, y => 0.0], [a => 1.0, b => 1.0], lcons = [1.0],
                           ucons = [1.0], grad = true, hess = true)
sol = solve(prob, IPNewton())
@test sol.minimum < 1.0
prob.f.cons(out, sol.minimizer, [1.0, 1.0])
@test out ≈ [1.0]
sol = solve(prob, Ipopt.Optimizer())
@test sol.minimum < 1.0
prob.f.cons(out, sol.minimizer, [1.0, 1.0])
@test out ≈ [1.0]
sol = solve(prob, AmplNLWriter.Optimizer(Ipopt_jll.amplexe))
@test sol.minimum < 1.0
prob.f.cons(out, sol.minimizer, [1.0, 1.0])
@test out ≈ [1.0]

rosenbrock(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
x0 = zeros(2)
_p = [1.0, 100.0]

f = OptimizationFunction(rosenbrock, Optimization.AutoModelingToolkit())
prob = OptimizationProblem(f, x0, _p)
sol = solve(prob, Newton())

# issue #819
@testset "Combined system name collisions" begin
    sys2 = OptimizationSystem(loss, [x, y], [a, b], name = :sys1)
    @test_throws ArgumentError OptimizationSystem(loss2, [z], [β], systems = [sys1, sys2])
end

# observed variable handling
@variables OBS
@named sys2 = OptimizationSystem(loss, [x, y], [a, b]; observed = [OBS ~ x + y])
OBS2 = OBS
@test isequal(OBS2, @nonamespace sys2.OBS)
@unpack OBS = sys2
@test isequal(OBS2, OBS)

# nested constraints
@testset "nested systems" begin
    @variables x y
    o1 = (x - 1)^2
    o2 = (y - 1 / 2)^2
    c1 = [
        x ~ 1,
    ]
    c2 = [
        y ~ 1,
    ]
    sys1 = OptimizationSystem(o1, [x], [], name = :sys1, constraints = c1)
    sys2 = OptimizationSystem(o2, [y], [], name = :sys2, constraints = c2)
    sys = OptimizationSystem(0, [], []; name = :sys, systems = [sys1, sys2],
                             constraints = [sys1.x + sys2.y ~ 2], checks = false)
    prob = OptimizationProblem(sys, [0.0, 0.0])

    @test isequal(constraints(sys), vcat(sys1.x + sys2.y ~ 2, sys1.x ~ 1, sys2.y ~ 1))
    @test isequal(equations(sys), (sys1.x - 1)^2 + (sys2.y - 1 / 2)^2)
    @test isequal(states(sys), [sys1.x, sys2.y])
end

@testset "time dependent var" begin
    @parameters t
    @variables x(t) y
    @parameters a b
    loss = (a - x)^2 + b * (y - x^2)^2
    sys1 = OptimizationSystem(loss, [x, y], [a, b], name = :sys1)

    cons2 = [x^2 + y^2 ~ 0, y * sin(x) - x ~ 0]
    sys2 = OptimizationSystem(loss, [x, y], [a, b], name = :sys2, constraints = cons2)

    @variables z
    @parameters β
    loss2 = sys1.x - sys2.y + z * β
    combinedsys = OptimizationSystem(loss2, [z], [β], systems = [sys1, sys2],
                                     name = :combinedsys)

    u0 = [sys1.x => 1.0
          sys1.y => 2.0
          sys2.x => 3.0
          sys2.y => 4.0
          z => 5.0]
    p = [sys1.a => 6.0
         sys1.b => 7.0
         sys2.a => 8.0
         sys2.b => 9.0
         β => 10.0]

    prob = OptimizationProblem(combinedsys, u0, p, grad = true)
    @test prob.f.sys === combinedsys
    sol = solve(prob, NelderMead())
    @test sol.minimum < -1e5

    prob2 = remake(prob, u0 = sol.minimizer)
    sol = solve(prob, BFGS(initial_stepnorm = 0.0001), allow_f_increases = true)
    @test sol.minimum < -1e8
    sol = solve(prob2, BFGS(initial_stepnorm = 0.0001), allow_f_increases = true)
    @test sol.minimum < -1e9

    #inequality constraint, the bounds for constraints lcons !== ucons
    prob = OptimizationProblem(sys2, [x => 0.0, y => 0.0], [a => 1.0, b => 100.0],
                               lcons = [-1.0, -1.0], ucons = [500.0, 500.0], grad = true,
                               hess = true)
    @test prob.f.sys === sys2
    sol = solve(prob, IPNewton(), allow_f_increases = true)
    @test sol.minimum < 1.0
    sol = solve(prob, Ipopt.Optimizer())
    @test sol.minimum < 1.0
end

@variables x
o1 = (x - 1)^2
c1 = [
    x ~ 1,
]
testdict = Dict(["test" => 1])
sys1 = OptimizationSystem(o1, [x], [], name = :sys1, constraints = c1,
                          metadata = testdict)
@test get_metadata(sys1) == testdict
