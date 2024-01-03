using ModelingToolkit, SparseArrays, Test, Optimization, OptimizationOptimJL,
    OptimizationMOI, Ipopt, AmplNLWriter, Ipopt_jll
using ModelingToolkit: get_metadata

@testset "basic" begin
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

    prob = OptimizationProblem(combinedsys, u0, p, grad = true, hess = true, cons_j = true,
        cons_h = true)
    @test prob.f.sys === combinedsys
    sol = solve(prob, Ipopt.Optimizer(); print_level = 0)
    @test sol.minimum < -1e5
end

@testset "inequality constraint" begin
    @variables x y
    @parameters a b
    loss = (a - x)^2 + b * (y - x^2)^2
    cons = [
        x^2 + y^2 ≲ 1.0,
    ]
    @named sys = OptimizationSystem(loss, [x, y], [a, b], constraints = cons)

    prob = OptimizationProblem(sys, [x => 0.0, y => 0.0], [a => 1.0, b => 1.0],
        grad = true, hess = true, cons_j = true, cons_h = true)
    @test prob.f.sys === sys
    sol = solve(prob, IPNewton())
    @test sol.minimum < 1.0
    sol = solve(prob, Ipopt.Optimizer(); print_level = 0)
    @test sol.minimum < 1.0

    prob = OptimizationProblem(sys, [x => 0.0, y => 0.0], [a => 1.0, b => 1.0],
        grad = false, hess = false, cons_j = false, cons_h = false)
    sol = solve(prob, AmplNLWriter.Optimizer(Ipopt_jll.amplexe))
    @test_skip sol.minimum < 1.0
end

@testset "equality constraint" begin
    @variables x y z
    @parameters a b
    loss = (a - x)^2 + b * z^2
    cons = [1.0 ~ x^2 + y^2
        z ~ y - x^2
        z^2 + y^2 ≲ 1.0]
    @named sys = OptimizationSystem(loss, [x, y, z], [a, b], constraints = cons)
    sys = structural_simplify(sys)
    prob = OptimizationProblem(sys, [x => 0.0, y => 0.0, z => 0.0], [a => 1.0, b => 1.0],
        grad = true, hess = true, cons_j = true, cons_h = true)
    sol = solve(prob, IPNewton())
    @test sol.minimum < 1.0
    @test sol.u≈[0.808, -0.064] atol=1e-3
    @test sol[x]^2 + sol[y]^2 ≈ 1.0
    sol = solve(prob, Ipopt.Optimizer(); print_level = 0)
    @test sol.minimum < 1.0
    @test sol.u≈[0.808, -0.064] atol=1e-3
    @test sol[x]^2 + sol[y]^2 ≈ 1.0

    prob = OptimizationProblem(sys, [x => 0.0, y => 0.0, z => 0.0], [a => 1.0, b => 1.0],
        grad = false, hess = false, cons_j = false, cons_h = false)
    sol = solve(prob, AmplNLWriter.Optimizer(Ipopt_jll.amplexe))
    @test_skip sol.minimum < 1.0
    @test_skip sol.u≈[0.808, -0.064] atol=1e-3
    @test_skip sol[x]^2 + sol[y]^2 ≈ 1.0
end

@testset "rosenbrock" begin
    rosenbrock(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
    x0 = zeros(2)
    p = [1.0, 100.0]
    f = OptimizationFunction(rosenbrock, Optimization.AutoModelingToolkit())
    prob = OptimizationProblem(f, x0, p)
    sol = solve(prob, Newton())
    @test sol.u ≈ [1.0, 1.0]
end

# issue #819
@testset "Combined system name collisions" begin
    @variables x y
    @parameters a b
    loss = (a - x)^2 + b * (y - x^2)^2
    sys1 = OptimizationSystem(loss, [x, y], [a, b], name = :sys1)
    sys2 = OptimizationSystem(loss, [x, y], [a, b], name = :sys1)
    @variables z
    @parameters β
    loss2 = sys1.x - sys2.y + z * β
    @test_throws ArgumentError OptimizationSystem(loss2, [z], [β], systems = [sys1, sys2])
end

@testset "observed variable handling" begin
    @variables x y
    @parameters a b
    loss = (a - x)^2 + b * (y - x^2)^2
    @variables OBS
    @named sys2 = OptimizationSystem(loss, [x, y], [a, b]; observed = [OBS ~ x + y])
    OBS2 = OBS
    @test isequal(OBS2, @nonamespace sys2.OBS)
    @unpack OBS = sys2
    @test isequal(OBS2, OBS)
end

# nested constraints
@testset "nested systems" begin
    @variables x y
    @parameters a = 1
    o1 = (x - a)^2
    o2 = (y - 1 / 2)^2
    c1 = [
        x ~ 1,
    ]
    c2 = [
        y ~ 1,
    ]
    sys1 = OptimizationSystem(o1, [x], [a], name = :sys1, constraints = c1)
    sys2 = OptimizationSystem(o2, [y], [], name = :sys2, constraints = c2)
    sys = OptimizationSystem(0, [], []; name = :sys, systems = [sys1, sys2],
        constraints = [sys1.x + sys2.y ~ 2], checks = false)
    prob = OptimizationProblem(sys, [0.0, 0.0])
    @test isequal(constraints(sys), vcat(sys1.x + sys2.y ~ 2, sys1.x ~ 1, sys2.y ~ 1))
    @test isequal(equations(sys), (sys1.x - sys1.a)^2 + (sys2.y - 1 / 2)^2)
    @test isequal(states(sys), [sys1.x, sys2.y])

    prob_ = remake(prob, u0 = [1.0, 0.0], p = [2.0])
    @test isequal(prob_.u0, [1.0, 0.0])
    @test isequal(prob_.p, [2.0])

    prob_ = remake(prob, u0 = Dict(sys1.x => 1.0), p = Dict(sys1.a => 2.0))
    @test isequal(prob_.u0, [1.0, 0.0])
    @test isequal(prob_.p, [2.0])
end

@testset "nested systems" begin
    @variables x1 x2 x3 x4
    @named sys1 = OptimizationSystem(x1, [x1], [])
    @named sys2 = OptimizationSystem(x2, [x2], [], systems = [sys1])
    @named sys3 = OptimizationSystem(x3, [x3], [], systems = [sys2])
    @named sys4 = OptimizationSystem(x4, [x4], [], systems = [sys3])

    @test isequal(equations(sys4), sys3.sys2.sys1.x1 + sys3.sys2.x2 + sys3.x3 + x4)
end

@testset "time dependent var" begin
    @parameters t
    @variables x(t) y
    @parameters a b
    loss = (a - x)^2 + b * (y - x^2)^2
    sys1 = OptimizationSystem(loss, [x, y], [a, b], name = :sys1)

    cons = [
        x^2 + y^2 ≲ 1.0,
    ]
    sys2 = OptimizationSystem(loss, [x, y], [a, b], name = :sys2, constraints = cons)

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

    prob = OptimizationProblem(combinedsys, u0, p, grad = true, hess = true, cons_j = true,
        cons_h = true)
    @test prob.f.sys === combinedsys
    @test_broken SciMLBase.successful_retcode(solve(prob,
        Ipopt.Optimizer();
        print_level = 0))
    #=
     @test sol.minimum < -1e5

     prob = OptimizationProblem(sys2, [x => 0.0, y => 0.0], [a => 1.0, b => 100.0],
         grad = true, hess = true, cons_j = true, cons_h = true)
     @test prob.f.sys === sys2
     sol = solve(prob, IPNewton())
     @test sol.minimum < 1.0
     sol = solve(prob, Ipopt.Optimizer(); print_level = 0)
     @test sol.minimum < 1.0
     =#
end

@testset "metadata" begin
    @variables x
    o1 = (x - 1)^2
    c1 = [
        x ~ 1,
    ]
    testdict = Dict(["test" => 1])
    sys1 = OptimizationSystem(o1, [x], [], name = :sys1, constraints = c1,
        metadata = testdict)
    @test get_metadata(sys1) == testdict
end

@testset "non-convex problem with inequalities" begin
    @variables x[1:2] [bounds = (0.0, Inf)]
    @named sys = OptimizationSystem(x[1] + x[2], [x...], [];
        constraints = [
            1.0 ≲ x[1]^2 + x[2]^2,
            x[1]^2 + x[2]^2 ≲ 2.0,
        ])

    prob = OptimizationProblem(sys, [x[1] => 2.0, x[2] => 0.0], [], grad = true,
        hess = true, cons_j = true, cons_h = true)
    sol = Optimization.solve(prob, Ipopt.Optimizer(); print_level = 0)
    @test sol.u ≈ [1, 0]
    @test prob.lb == [0.0, 0.0]
    @test prob.ub == [Inf, Inf]
end

@testset "parameter bounds" begin
    @parameters c = 0.0
    @variables x y [bounds = (c, Inf)]
    @parameters a b
    loss = (a - x)^2 + b * (y - x^2)^2
    @named sys = OptimizationSystem(loss, [x, y], [a, b, c])
    prob = OptimizationProblem(sys, [x => 0.0, y => 0.0], [a => 1.0, b => 100.0])
    @test prob.lb == [-Inf, 0.0]
    @test prob.ub == [Inf, Inf]
end

@testset "modelingtoolkitize" begin
    @variables x₁ x₂
    @parameters α₁ α₂
    loss = (α₁ - x₁)^2 + α₂ * (x₂ - x₁^2)^2
    cons = [
        x₁^2 + x₂^2 ≲ 1.0,
    ]
    sys1 = OptimizationSystem(loss, [x₁, x₂], [α₁, α₂], name = :sys1, constraints = cons)

    prob1 = OptimizationProblem(sys1, [x₁ => 0.0, x₂ => 0.0], [α₁ => 1.0, α₂ => 100.0],
        grad = true, hess = true, cons_j = true, cons_h = true)

    sys2 = modelingtoolkitize(prob1)
    prob2 = OptimizationProblem(sys2, [x₁ => 0.0, x₂ => 0.0], [α₁ => 1.0, α₂ => 100.0],
        grad = true, hess = true, cons_j = true, cons_h = true)

    sol1 = Optimization.solve(prob1, Ipopt.Optimizer())
    sol2 = Optimization.solve(prob2, Ipopt.Optimizer())

    @test sol1.u ≈ sol2.u
end

@testset "#2323 keep symbolic expressions and xor condition on constraint bounds" begin
    @variables x y
    @parameters a b
    loss = (a - x)^2 + b * (y - x^2)^2
    @named sys = OptimizationSystem(loss, [x, y], [a, b], constraints = [x^2 + y^2 ≲ 0.0])
    @test_throws ArgumentError OptimizationProblem(sys,
        [x => 0.0, y => 0.0],
        [a => 1.0, b => 100.0],
        lcons = [0.0])
    @test_throws ArgumentError OptimizationProblem(sys,
        [x => 0.0, y => 0.0],
        [a => 1.0, b => 100.0],
        ucons = [0.0])

    prob = OptimizationProblem(sys, [x => 0.0, y => 0.0], [a => 1.0, b => 100.0])
    @test prob.f.expr isa Symbolics.Symbolic
    @test all(prob.f.cons_expr[i].lhs isa Symbolics.Symbolic
              for i in 1:length(prob.f.cons_expr))
end
