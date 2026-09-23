using ModelingToolkitBase, Test
using ModelingToolkitBase: costs, cost, accepts_array_equations
using SciMLBase
using OptimizationEvolutionary
using Symbolics
using LinearAlgebra
using StableRNGs

function moo_system()
    @variables x y
    @parameters a
    costs = [(x - a)^2, (y - 1)^2]
    sys = complete(OptimizationSystem(costs, [x, y], [a]; name = :moo))
    return sys, x, y, a
end

@testset "`costs` returns the unconsolidated objective vector" begin
    sys, x, y, a = moo_system()
    objs = costs(sys)
    @test length(objs) == 2
    @test isequal(objs[1], (x - a)^2)
    @test isequal(objs[2], (y - 1)^2)
    # `cost` is the same vector folded through `consolidate`
    @test isequal(cost(sys), (x - a)^2 + (y - 1)^2)
end

@testset "MultiObjectiveOptimizationFunction constructs and evaluates" begin
    sys, x, y, a = moo_system()

    f = MultiObjectiveOptimizationFunction{false}(sys)
    @test f isa MultiObjectiveOptimizationFunction{false}
    @test f.f([1.5, 2.0], [3.0]) ≈ [(1.5 - 3.0)^2, (2.0 - 1)^2]

    fiip = MultiObjectiveOptimizationFunction{true}(sys)
    @test fiip isa MultiObjectiveOptimizationFunction{true}
    out = zeros(2)
    fiip.f(out, [1.5, 2.0], [3.0])
    @test out ≈ [(1.5 - 3.0)^2, (2.0 - 1)^2]

    # the scalar `OptimizationFunction` path still consolidates
    fscalar = OptimizationFunction{false}(sys)
    @test fscalar.f([1.5, 2.0], [3.0]) ≈ (1.5 - 3.0)^2 + (2.0 - 1)^2
end

@testset "objective jacobian and hessians" begin
    sys, x, y, a = moo_system()
    f = MultiObjectiveOptimizationFunction{false}(sys; jac = true, hess = true)

    u = [1.5, 2.0]
    p = [3.0]
    @test f.jac(u, p) ≈ [2 * (u[1] - p[1]) 0.0; 0.0 2 * (u[2] - 1.0)]

    # one hessian function per objective, as `OptimizationBase` consumes `hess`
    @test length(f.hess) == 2
    @test f.hess[1](u, p) ≈ [2.0 0.0; 0.0 0.0]
    @test f.hess[2](u, p) ≈ [0.0 0.0; 0.0 2.0]

    fiip = MultiObjectiveOptimizationFunction{true}(sys; jac = true, hess = true)
    J = zeros(2, 2)
    fiip.jac(J, u, p)
    @test J ≈ f.jac(u, p)
    H = zeros(2, 2)
    fiip.hess[2](H, u, p)
    @test H ≈ [0.0 0.0; 0.0 2.0]

    fsp = MultiObjectiveOptimizationFunction{false}(sys; hess = true, sparse = true)
    @test length(fsp.hess_prototype) == 2
    @test fsp.hess[1](u, p) ≈ [2.0 0.0; 0.0 0.0]

    fexpr = MultiObjectiveOptimizationFunction{false}(
        sys; jac = true, hess = true, expression = Val{true}
    )
    @test fexpr isa Expr
end

@testset "constraints are shared with the scalar path" begin
    @variables x y
    sys = complete(
        OptimizationSystem(
            [(x - 1)^2, (y - 2)^2], [x, y], [];
            constraints = [x + y ≲ 3.0], name = :cmoo
        )
    )
    f = MultiObjectiveOptimizationFunction{false}(sys; cons_j = true)
    @test f.cons !== nothing
    @test f.cons([1.0, 1.0], nothing) ≈ [1.0 + 1.0 - 3.0]
    @test f.cons_j([1.0, 1.0], nothing) ≈ [1.0 1.0]

    prob = OptimizationProblem(sys, [x => 0.0, y => 0.0]; multiobjective = true)
    @test prob.lcons == [-Inf]
    @test prob.ucons == [0.0]
end

@testset "`OptimizationProblem` routes through `multiobjective`" begin
    sys, x, y, a = moo_system()

    prob = OptimizationProblem(
        sys, [x => 0.0, y => 0.0, a => 3.0]; multiobjective = true
    )
    @test prob.f isa MultiObjectiveOptimizationFunction
    @test prob.f.f([1.5, 2.0], prob.p) ≈ [(1.5 - 3.0)^2, (2.0 - 1)^2]

    # bounds machinery still applies
    prob_b = OptimizationProblem(
        sys, [x => 0.0, y => 0.0, a => 3.0];
        multiobjective = true, lb = [-1.0, -1.0], ub = [1.0, 1.0]
    )
    @test prob_b.lb == [-1.0, -1.0]
    @test prob_b.ub == [1.0, 1.0]

    # the default path is unchanged
    prob_scalar = OptimizationProblem(sys, [x => 0.0, y => 0.0, a => 3.0])
    @test prob_scalar.f isa OptimizationFunction
    @test prob_scalar.f.f([1.5, 2.0], prob_scalar.p) ≈ (1.5 - 3.0)^2 + (2.0 - 1)^2

    # `weights` scalarizes, so it cannot combine with `multiobjective`
    @test_throws ArgumentError OptimizationProblem(
        sys, [x => 0.0, y => 0.0, a => 3.0];
        multiobjective = true, weights = [1.0, 1.0]
    )
    @test_throws ArgumentError MultiObjectiveOptimizationFunction{false}(
        sys; weights = [1.0, 1.0]
    )
end

@testset "array unknowns and weighted sums" begin
    @variables w[1:2]
    @parameters a
    sys = complete(
        System(
            Equation[], [w], [a];
            costs = [sum(abs2, w .- a), sum(abs2, w .+ a)], name = :array_moo
        )
    )
    @test length(unknowns(sys)) == 1
    op = [w => [1.0, 2.0], a => 3.0]
    multi = OptimizationProblem(sys, op; multiobjective = true, jac = true)
    @test multi.u0 == [1.0, 2.0]
    @test multi.f.f(multi.u0, multi.p) ≈ [5.0, 41.0]
    # rows are costs, columns are flattened unknowns: ∂/∂w of each cost at w=[1,2], a=3
    @test multi.f.jac(multi.u0, multi.p) ≈ [-4.0 -2.0; 8.0 10.0]

    weights = [0.25, 0.75]
    weighted = OptimizationProblem(sys, op; weights)
    @test weighted.f(weighted.u0, weighted.p) ≈ dot(weights, multi.f.f(multi.u0, multi.p))
end

@testset "multi-objective NSGA-II consumes generated costs" begin
    # Solver-independent check that NSGA-II evaluates the generated multi-objective
    # function. No Pareto-front claim: a seed sweep of the analytic front check fails
    # for most streams, and the unbounded path seeds the population from u0.
    @variables x
    sys = complete(
        System(
            Equation[], [x], []; costs = [x^2, (x - 2)^2], name = :pareto
        )
    )
    # Off-front u0 + bounds so the initial population is random, not copies of a
    # Pareto point. Requires OptimizationEvolutionary ≥ 0.4.13 (bounded MOO path).
    prob = OptimizationProblem(
        sys, [x => 5.0]; multiobjective = true, lb = [-10.0], ub = [10.0]
    )
    sol = solve(
        prob, OptimizationEvolutionary.NSGA2(); maxiters = 40, rng = StableRNG(1234)
    )

    @test !isempty(sol.u)
    for u in sol.u
        objectives = prob.f.f(u, prob.p)
        @test objectives ≈ [u[1]^2, (u[1] - 2)^2]
    end
end

@testset "equations are still rejected" begin
    # `MultiObjectiveOptimizationFunction` vectorizes `costs`, not `equations`:
    # scalar and array equations alike are rejected by `check_no_equations`.
    @test !accepts_array_equations(MultiObjectiveOptimizationFunction)

    @variables x u[1:2]
    @named esys = System([x ~ 2.0], [x], []; costs = [x^2])
    @test_throws ["cannot be used"] MultiObjectiveOptimizationFunction{false}(
        complete(esys)
    )

    @named asys = System([u ~ [1.0, 2.0]], [u...], []; costs = [sum(u)])
    @test_throws ["cannot be used"] MultiObjectiveOptimizationFunction{false}(
        complete(asys)
    )
    @test_throws ["cannot be used"] OptimizationProblem(
        complete(asys), [u[1] => 0.0, u[2] => 0.0]; multiobjective = true
    )
end
