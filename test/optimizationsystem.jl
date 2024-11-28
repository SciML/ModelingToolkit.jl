using ModelingToolkit, SparseArrays, Test, Optimization, OptimizationOptimJL,
      OptimizationMOI, Ipopt, AmplNLWriter, Ipopt_jll, SymbolicIndexingInterface
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
    combinedsys = complete(OptimizationSystem(loss2, [z], [β], systems = [sys1, sys2],
        name = :combinedsys))

    equations(combinedsys)
    unknowns(combinedsys)
    parameters(combinedsys)

    calculate_gradient(combinedsys)
    calculate_hessian(combinedsys)
    generate_function(combinedsys)
    generate_gradient(combinedsys)
    generate_hessian(combinedsys)
    hess_sparsity = ModelingToolkit.hessian_sparsity(sys1)
    sparse_prob = OptimizationProblem(complete(sys1),
        [x, y],
        [a => 0.0, b => 0.0],
        grad = true,
        sparse = true)
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
    @test sol.objective < -1e5
end

@testset "inequality constraint" begin
    @variables x y
    @parameters a b
    loss = (a - x)^2 + b * (y - x^2)^2
    cons = [
        x^2 + y^2 ≲ 1.0
    ]
    @named sys = OptimizationSystem(loss, [x, y], [a, b], constraints = cons)
    sys = complete(sys)
    prob = OptimizationProblem(sys, [x => 0.0, y => 0.0], [a => 1.0, b => 1.0],
        grad = true, hess = true, cons_j = true, cons_h = true)
    @test prob.f.sys === sys
    sol = solve(prob, IPNewton())
    @test sol.objective < 1.0
    sol = solve(prob, Ipopt.Optimizer(); print_level = 0)
    @test sol.objective < 1.0

    prob = OptimizationProblem(sys, [x => 0.0, y => 0.0], [a => 1.0, b => 1.0],
        grad = false, hess = false, cons_j = false, cons_h = false)
    sol = solve(prob, AmplNLWriter.Optimizer(Ipopt_jll.amplexe))
    @test_skip sol.objective < 1.0
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
    @test sol.objective < 1.0
    @test sol.u≈[0.808, -0.064] atol=1e-3
    @test sol[x]^2 + sol[y]^2 ≈ 1.0
    sol = solve(prob, Ipopt.Optimizer(); print_level = 0)
    @test sol.objective < 1.0
    @test sol.u≈[0.808, -0.064] atol=1e-3
    @test sol[x]^2 + sol[y]^2 ≈ 1.0

    prob = OptimizationProblem(sys, [x => 0.0, y => 0.0, z => 0.0], [a => 1.0, b => 1.0],
        grad = false, hess = false, cons_j = false, cons_h = false)
    sol = solve(prob, AmplNLWriter.Optimizer(Ipopt_jll.amplexe))
    @test_skip sol.objective < 1.0
    @test_skip sol.u≈[0.808, -0.064] atol=1e-3
    @test_skip sol[x]^2 + sol[y]^2 ≈ 1.0
end

@testset "rosenbrock" begin
    rosenbrock(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
    x0 = zeros(2)
    p = [1.0, 100.0]
    f = OptimizationFunction(rosenbrock, Optimization.AutoSymbolics())
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
        x ~ 1
    ]
    c2 = [
        y ~ 1
    ]
    sys1 = OptimizationSystem(o1, [x], [a], name = :sys1, constraints = c1)
    sys2 = OptimizationSystem(o2, [y], [], name = :sys2, constraints = c2)
    sys = complete(OptimizationSystem(0, [], []; name = :sys, systems = [sys1, sys2],
        constraints = [sys1.x + sys2.y ~ 2], checks = false))
    prob = OptimizationProblem(sys, [0.0, 0.0])
    @test isequal(constraints(sys), vcat(sys1.x + sys2.y ~ 2, sys1.x ~ 1, sys2.y ~ 1))
    @test isequal(equations(sys), (sys1.x - sys1.a)^2 + (sys2.y - 1 / 2)^2)
    @test isequal(unknowns(sys), [sys1.x, sys2.y])

    prob_ = remake(prob, u0 = [1.0, 0.0], p = [2.0])
    @test isequal(prob_.u0, [1.0, 0.0])
    @test isequal(prob_.p, [2.0])

    prob_ = remake(prob, u0 = Dict(sys1.x => 1.0), p = Dict(sys1.a => 2.0))
    @test isequal(prob_.u0, [1.0, 0.0])
    @test isequal((prob_.p...,)[1], [2.0])
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
    @independent_variables t
    @variables x(t) y
    @parameters a b
    loss = (a - x)^2 + b * (y - x^2)^2
    sys1 = OptimizationSystem(loss, [x, y], [a, b], name = :sys1)

    cons = [
        x^2 + y^2 ≲ 1.0
    ]
    sys2 = OptimizationSystem(loss, [x, y], [a, b], name = :sys2, constraints = cons)

    @variables z
    @parameters β
    loss2 = sys1.x - sys2.y + z * β
    combinedsys = complete(OptimizationSystem(loss2, [z], [β], systems = [sys1, sys2],
        name = :combinedsys))

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
     @test sol.objective < -1e5

     prob = OptimizationProblem(sys2, [x => 0.0, y => 0.0], [a => 1.0, b => 100.0],
         grad = true, hess = true, cons_j = true, cons_h = true)
     @test prob.f.sys === sys2
     sol = solve(prob, IPNewton())
     @test sol.objective < 1.0
     sol = solve(prob, Ipopt.Optimizer(); print_level = 0)
     @test sol.objective < 1.0
     =#
end

@testset "metadata" begin
    @variables x
    o1 = (x - 1)^2
    c1 = [
        x ~ 1
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
            x[1]^2 + x[2]^2 ≲ 2.0
        ])

    prob = OptimizationProblem(complete(sys), [x[1] => 2.0, x[2] => 0.0], [], grad = true,
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
    prob = OptimizationProblem(complete(sys), [x => 0.0, y => 0.0], [a => 1.0, b => 100.0])
    @test prob.lb == [-Inf, 0.0]
    @test prob.ub == [Inf, Inf]
end

@testset "modelingtoolkitize" begin
    @variables x₁ x₂
    @parameters α₁ α₂
    loss = (α₁ - x₁)^2 + α₂ * (x₂ - x₁^2)^2
    cons = [
        x₁^2 + x₂^2 ≲ 1.0
    ]
    sys1 = complete(OptimizationSystem(loss,
        [x₁, x₂],
        [α₁, α₂],
        name = :sys1,
        constraints = cons))

    prob1 = OptimizationProblem(sys1, [x₁ => 0.0, x₂ => 0.0], [α₁ => 1.0, α₂ => 100.0],
        grad = true, hess = true, cons_j = true, cons_h = true)

    sys2 = complete(modelingtoolkitize(prob1))
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
    sys = complete(sys)
    @test_throws Any OptimizationProblem(sys,
        [x => 0.0, y => 0.0],
        [a => 1.0, b => 100.0],
        lcons = [0.0])
    @test_throws Any OptimizationProblem(sys,
        [x => 0.0, y => 0.0],
        [a => 1.0, b => 100.0],
        ucons = [0.0])

    prob = OptimizationProblem(sys, [x => 0.0, y => 0.0], [a => 1.0, b => 100.0])
    @test prob.f.expr isa Symbolics.Symbolic
    @test all(prob.f.cons_expr[i].lhs isa Symbolics.Symbolic
    for i in 1:length(prob.f.cons_expr))
end

@testset "Derivatives, iip and oop" begin
    @variables x y
    @parameters a b
    loss = (a - x)^2 + b * (y - x^2)^2
    cons2 = [x^2 + y^2 ~ 0, y * sin(x) - x ~ 0]
    sys = complete(OptimizationSystem(
        loss, [x, y], [a, b], name = :sys2, constraints = cons2))
    prob = OptimizationProblem(sys, [x => 0.0, y => 0.0], [a => 1.0, b => 100.0],
        grad = true, hess = true, cons_j = true, cons_h = true)

    G1 = Array{Float64}(undef, 2)
    H1 = Array{Float64}(undef, 2, 2)
    J = Array{Float64}(undef, 2, 2)
    H3 = [Array{Float64}(undef, 2, 2), Array{Float64}(undef, 2, 2)]

    prob.f.grad(G1, [1.0, 1.0], [1.0, 100.0])
    @test prob.f.grad([1.0, 1.0], [1.0, 100.0]) == G1

    prob.f.hess(H1, [1.0, 1.0], [1.0, 100.0])
    @test prob.f.hess([1.0, 1.0], [1.0, 100.0]) == H1

    prob.f.cons_j(J, [1.0, 1.0], [1.0, 100.0])
    @test prob.f.cons_j([1.0, 1.0], [1.0, 100.0]) == J

    prob.f.cons_h(H3, [1.0, 1.0], [1.0, 100.0])
    @test prob.f.cons_h([1.0, 1.0], [1.0, 100.0]) == H3
end

@testset "Passing `nothing` to `u0`" begin
    @variables x = 1.0
    @mtkbuild sys = OptimizationSystem((x - 3)^2, [x], [])
    prob = @test_nowarn OptimizationProblem(sys, nothing)
    @test_nowarn solve(prob, NelderMead())
end

@testset "Bounded unknowns are irreducible" begin
    @variables x
    @variables y [bounds = (-Inf, Inf)]
    @variables z [bounds = (1.0, 2.0)]
    obj = x^2 + y^2 + z^2
    cons = [y ~ 2x
            z ~ 2y]
    @mtkbuild sys = OptimizationSystem(obj, [x, y, z], []; constraints = cons)
    @test is_variable(sys, z)
    @test !is_variable(sys, y)

    @variables x[1:3] [bounds = ([-Inf, -1.0, -2.0], [Inf, 1.0, 2.0])]
    obj = x[1]^2 + x[2]^2 + x[3]^2
    cons = [x[2] ~ 2x[1] + 3, x[3] ~ x[1] + x[2]]
    @mtkbuild sys = OptimizationSystem(obj, [x], []; constraints = cons)
    @test length(unknowns(sys)) == 2
    @test !is_variable(sys, x[1])
    @test is_variable(sys, x[2])
    @test is_variable(sys, x[3])
end

@testset "Constraints work with nonnumeric parameters" begin
    @variables x
    @parameters p f(::Real)
    @mtkbuild sys = OptimizationSystem(
        x^2 + f(x) * p, [x], [f, p]; constraints = [2.0 ≲ f(x) + p])
    prob = OptimizationProblem(sys, [x => 1.0], [p => 1.0, f => (x -> 2x)])
    @test abs(prob.f.cons(prob.u0, prob.p)[1]) ≈ 1.0
end

@testset "Variable discovery" begin
    @variables x1 x2
    @parameters p1 p2
    @named sys1 = OptimizationSystem(x1^2; constraints = [p1 * x1 ≲ 2.0])
    @named sys2 = OptimizationSystem(x2^2; constraints = [p2 * x2 ≲ 2.0], systems = [sys1])
    @test isequal(only(unknowns(sys1)), x1)
    @test isequal(only(parameters(sys1)), p1)
    @test all(y -> any(x -> isequal(x, y), unknowns(sys2)), [x2, sys1.x1])
    @test all(y -> any(x -> isequal(x, y), parameters(sys2)), [p2, sys1.p1])
end
