using ModelingToolkit, StaticArrays, LinearAlgebra
using DiffEqBase, SparseArrays
using Test
using NonlinearSolve
using ForwardDiff
using SymbolicIndexingInterface
using ModelingToolkit: value
using ModelingToolkit: get_default_or_guess, MTKParameters

canonequal(a, b) = isequal(simplify(a), simplify(b))

# Define some variables
@parameters σ ρ β
@constants h = 1
@variables x y z

function test_nlsys_inference(name, sys, vs, ps)
    @testset "NonlinearSystem construction: $name" begin
        @test Set(unknowns(sys)) == Set(value.(vs))
        @test Set(parameters(sys)) == Set(value.(ps))
    end
end

# Define a nonlinear system
eqs = [0 ~ σ * (y - x) * h,
    0 ~ x * (ρ - z) - y,
    0 ~ x * y - β * z]
@named ns = System(eqs, [x, y, z], [σ, ρ, β, h], defaults = Dict(x => 2))
ns2 = eval(toexpr(ns))
@test issetequal(equations(ns), equations(ns2))
@test issetequal(unknowns(ns), unknowns(ns2))
@test issetequal(parameters(ns), parameters(ns2))
test_nlsys_inference("standard", ns, (x, y, z), (σ, ρ, β, h))
@test begin
    f = generate_rhs(ns, expression = Val{false})[2]
    du = [0.0, 0.0, 0.0]
    f(du, [1, 2, 3], [1, 2, 3, 1])
    du ≈ [1, -3, -7]
end

# Now nonlinear system with only variables
@variables x y z
@parameters σ ρ β

# Define a nonlinear system
eqs = [0 ~ σ * (y - x),
    y ~ x * (ρ - z),
    β * z ~ x * y]
@named ns = System(eqs, [x, y, z], [σ, ρ, β])
jac = calculate_jacobian(ns)
@testset "nlsys jacobian" begin
    @test canonequal(jac[1, 1], σ * -1)
    @test canonequal(jac[1, 2], σ)
    @test canonequal(jac[1, 3], 0)
    @test canonequal(jac[2, 1], ρ - z)
    @test canonequal(jac[2, 2], -1)
    @test canonequal(jac[2, 3], x * -1)
    @test canonequal(jac[3, 1], y)
    @test canonequal(jac[3, 2], x)
    @test canonequal(jac[3, 3], -1 * β)
end

# Intermediate calculations
a = y - x
# Define a nonlinear system
eqs = [0 ~ σ * a * h,
    0 ~ x * (ρ - z) - y,
    0 ~ x * y - β * z]
@named ns = System(eqs, [x, y, z], [σ, ρ, β, h])
ns = complete(ns)
nlsys_func = generate_rhs(ns)
nf = NonlinearFunction(ns)
jac = calculate_jacobian(ns)

@test ModelingToolkit.jacobian_sparsity(ns).colptr == sparse(jac).colptr
@test ModelingToolkit.jacobian_sparsity(ns).rowval == sparse(jac).rowval

jac = generate_jacobian(ns)

sH = calculate_hessian(ns)
@test getfield.(ModelingToolkit.hessian_sparsity(ns), :colptr) ==
      getfield.(sparse.(sH), :colptr)
@test getfield.(ModelingToolkit.hessian_sparsity(ns), :rowval) ==
      getfield.(sparse.(sH), :rowval)

prob = NonlinearProblem(ns, [x => 1.0, y => 1.0, z => 1.0, σ => 1.0, ρ => 1.0, β => 1.0])
@test prob.f.sys === ns
sol = solve(prob, NewtonRaphson())
@test sol.u[1] ≈ sol.u[2]

prob = NonlinearProblem(
    ns, [x => 1.0, y => 1.0, z => 1.0, σ => 1.0, ρ => 1.0, β => 1.0], jac = true)
@test_nowarn solve(prob, NewtonRaphson())

@variables u F s a
eqs1 = [
    0 ~ σ * (y - x) * h + F,
    0 ~ x * (ρ - z) - u,
    0 ~ x * y - β * z,
    0 ~ x + y - z - u
]

lorenz = name -> System(eqs1, [x, y, z, u, F], [σ, ρ, β, h], name = name)
lorenz1 = lorenz(:lorenz1)
@test_throws ArgumentError NonlinearProblem(complete(lorenz1), zeros(4))
lorenz2 = lorenz(:lorenz2)
@named connected = System(
    [s ~ a + lorenz1.x
     lorenz2.y ~ s * h
     lorenz1.F ~ lorenz2.u
     lorenz2.F ~ lorenz1.u],
    [s, a], [h],
    systems = [lorenz1, lorenz2])
@test_nowarn alias_elimination(connected)

# system promotion
using OrdinaryDiffEq
@independent_variables t
D = Differential(t)
@named subsys = convert_system_indepvar(lorenz1, t)
@named sys = System([D(subsys.x) ~ subsys.x + subsys.x], t, systems = [subsys])
sys = mtkcompile(sys)
u0 = [subsys.x => 1, subsys.z => 2.0, subsys.y => 1.0]
prob = ODEProblem(sys, [u0; [subsys.σ => 1, subsys.ρ => 2, subsys.β => 3]], (0, 1.0))
sol = solve(prob, FBDF(), reltol = 1e-7, abstol = 1e-7)
@test sol[subsys.x] + sol[subsys.y] - sol[subsys.z]≈sol[subsys.u] atol=1e-7
@test_throws ArgumentError convert_system_indepvar(sys, t)

@parameters σ ρ β
@variables x y z

# Define a nonlinear system
eqs = [0 ~ σ * (y - x),
    0 ~ x * (ρ - z) - y,
    0 ~ x * y - β * z * h]
@named ns = System(eqs, [x, y, z], [σ, ρ, β, h])
np = NonlinearProblem(
    complete(ns), [x => 0, y => 0, z => 0, σ => 1, ρ => 2, β => 3], jac = true, sparse = true)
@test calculate_jacobian(ns, sparse = true) isa SparseMatrixCSC

# issue #819
@testset "Combined system name collisions" begin
    function makesys(name)
        @parameters a
        @variables x f

        System([0 ~ -a * x + f], [x, f], [a]; name)
    end

    function issue819()
        sys1 = makesys(:sys1)
        sys2 = makesys(:sys1)
        @test_throws ModelingToolkit.NonUniqueSubsystemsError System(
            [sys2.f ~ sys1.x, sys1.f ~ 0], [], [],
            systems = [sys1, sys2], name = :foo)
    end
    issue819()
end

# issue #1115
@testset "Extending a NonlinearSystem with no iv" begin
    @parameters a b
    @variables x y
    eqs1 = [
        0 ~ a * x
    ]
    eqs2 = [
        0 ~ b * y
    ]

    @named sys1 = System(eqs1, [x], [a])
    @named sys2 = System(eqs2, [y], [b])
    @named sys3 = extend(sys1, sys2)

    @test isequal(union(Set(parameters(sys1)), Set(parameters(sys2))),
        Set(parameters(sys3)))
    @test isequal(union(Set(unknowns(sys1)), Set(unknowns(sys2))), Set(unknowns(sys3)))
    @test isequal(union(Set(equations(sys1)), Set(equations(sys2))), Set(equations(sys3)))
end

# observed variable handling
@independent_variables t
@parameters τ
@variables x(t) RHS(t)
@named fol = System([0 ~ (1 - x * h) / τ], [x], [τ];
    observed = [RHS ~ (1 - x) / τ])
@test isequal(RHS, @nonamespace fol.RHS)
RHS2 = RHS
@unpack RHS = fol
@test isequal(RHS, RHS2)

# issue #1358
@independent_variables t
@variables v1(t) v2(t) i1(t) i2(t)
eq = [v1 ~ sin(2pi * t * h)
      v1 - v2 ~ i1
      v2 ~ i2
      i1 ~ i2]
@named sys = System(eq, t)
@test length(equations(mtkcompile(sys))) == 0

@testset "Remake" begin
    @parameters a=1.0 b=1.0 c=1.0
    @constants h = 1
    @variables x y z

    eqs = [0 ~ a * (y - x) * h,
        0 ~ x * (b - z) - y,
        0 ~ x * y - c * z]
    @named sys = System(eqs, [x, y, z], [a, b, c, h], defaults = Dict(x => 2.0))
    sys = complete(sys)
    prob = NonlinearProblem(sys, ones(length(unknowns(sys))))

    prob_ = remake(prob, u0 = [1.0, 2.0, 3.0], p = [a => 1.1, b => 1.2, c => 1.3])
    @test prob_.u0 == [1.0, 2.0, 3.0]
    initials = unknowns(sys) .=> [1.0, 2.0, 3.0]
    @test prob_.p == MTKParameters(sys, [a => 1.1, b => 1.2, c => 1.3, initials...])

    prob_ = remake(prob, u0 = Dict(y => 2.0), p = Dict(a => 2.0))
    @test prob_.u0 == [1.0, 2.0, 1.0]
    initials = [x => 1.0, y => 2.0, z => 1.0]
    @test prob_.p == MTKParameters(sys, [a => 2.0, b => 1.0, c => 1.0, initials...])
end

@testset "Observed function generation without parameters" begin
    @variables x y z

    eqs = [0 ~ x + sin(y),
        0 ~ z - cos(x),
        0 ~ x * y]
    @named ns = System(eqs, [x, y, z], [])
    ns = complete(ns)
    vs = [unknowns(ns); parameters(ns)]
    ss_mtk = mtkcompile(ns)
    prob = NonlinearProblem(ss_mtk, vs .=> 1.0)
    sol = solve(prob)
    @test_nowarn sol[unknowns(ns)]
end

# Issue#2625
@parameters p d
@variables X(t)
alg_eqs = [0 ~ p - d * X]

sys = @test_nowarn System(alg_eqs; name = :name)
@test isequal(only(unknowns(sys)), X)
@test all(isequal.(parameters(sys), [p, d]))

# Over-determined sys
@variables u1 u2
@parameters u3 u4
eqs = [u3 ~ u1 + u2, u4 ~ 2 * (u1 + u2), u3 + u4 ~ 3 * (u1 + u2)]
@named ns = System(eqs, [u1, u2], [u3, u4])
sys = mtkcompile(ns; fully_determined = false)
@test length(unknowns(sys)) == 1

# Conservative
@variables X(t)
alg_eqs = [1 ~ 2X]
@named ns = System(alg_eqs)
sys = mtkcompile(ns)
@test length(equations(sys)) == 0
sys = mtkcompile(ns; conservative = true)
@test length(equations(sys)) == 1

# https://github.com/SciML/ModelingToolkit.jl/issues/2858
@testset "Jacobian/Hessian with observed equations that depend on unknowns" begin
    @variables x y z
    @parameters σ ρ β
    eqs = [0 ~ σ * (y - x)
           0 ~ x * (ρ - z) - y
           0 ~ x * y - β * z]
    guesses = [x => 1.0, z => 0.0]
    ps = [σ => 10.0, ρ => 26.0, β => 8 / 3]
    @mtkcompile ns = System(eqs)

    @test isequal(calculate_jacobian(ns), [(-1 - z + ρ)*σ -x*σ
                                           2x*(-z + ρ) -β-(x^2)])
    # solve without analytical jacobian
    prob = NonlinearProblem(ns, [guesses; ps])
    sol = solve(prob, NewtonRaphson())
    @test sol.retcode == ReturnCode.Success

    # solve with analytical jacobian
    prob = NonlinearProblem(ns, [guesses; ps], jac = true)
    sol = solve(prob, NewtonRaphson())
    @test sol.retcode == ReturnCode.Success

    # system that contains a chain of observed variables when simplified
    @variables x y z
    eqs = [0 ~ x^2 + 2z + y, z ~ y, y ~ x] # analytical solution x = y = z = 0 or -3
    @mtkcompile ns = System(eqs) # solve for y with observed chain z -> y -> x
    @test isequal(expand.(calculate_jacobian(ns)), [-3 // 2 - x;;])
    @test isequal(calculate_hessian(ns), [[-1;;]])
    prob = NonlinearProblem(ns, unknowns(ns) .=> -4.0) # give guess < -3 to reach -3
    sol = solve(prob, NewtonRaphson())
    @test sol[x] ≈ sol[y] ≈ sol[z] ≈ -3
end

@testset "Passing `nothing` to `u0`" begin
    @variables x = 1
    @mtkcompile sys = System([0 ~ x^2 - x^3 + 3])
    prob = @test_nowarn NonlinearProblem(sys, nothing)
    @test_nowarn solve(prob)
end

@testset "System of linear equations with vector variable" begin
    # 1st example in https://en.wikipedia.org/w/index.php?title=System_of_linear_equations&oldid=1247697953
    @variables x[1:3]
    A = [3 2 -1
         2 -2 4
         -1 1/2 -1]
    b = [1, -2, 0]
    @named sys = System(A * x ~ b, [x], [])
    sys = mtkcompile(sys)
    prob = NonlinearProblem(sys, unknowns(sys) .=> 0.0)
    sol = solve(prob)
    @test all(sol[x] .≈ A \ b)
end

@testset "resid_prototype when system has no unknowns and an equation" begin
    @variables x
    @parameters p
    @named sys = System([x ~ 1, x^2 - p ~ 0])
    for sys in [
        mtkcompile(sys, fully_determined = false),
        mtkcompile(sys, fully_determined = false, split = false)
    ]
        @test length(equations(sys)) == 1
        @test length(unknowns(sys)) == 0
        T = typeof(ForwardDiff.Dual(1.0))
        prob = NonlinearProblem(sys, [p => ForwardDiff.Dual(1.0)]; check_length = false)
        @test prob.f(Float64[], prob.p) isa Vector{T}
        @test prob.f.resid_prototype isa Vector{T}
        @test_nowarn solve(prob)
    end
end

@testset "IntervalNonlinearProblem" begin
    @variables x
    @parameters p
    @named nlsys = System([0 ~ x * x - p])

    for sys in [complete(nlsys), complete(nlsys; split = false)]
        prob = IntervalNonlinearProblem(sys, (0.0, 2.0), [p => 1.0])
        sol = @test_nowarn solve(prob, ITP())
        @test SciMLBase.successful_retcode(sol)
        @test_nowarn IntervalNonlinearProblem(
            sys, (0.0, 2.0), [p => 1.0]; expression = Val{true})
    end

    @variables y
    @mtkcompile sys = System([0 ~ x * x - p * x + p, 0 ~ x * y + p])
    @test_throws ["single unknown"] IntervalNonlinearProblem(sys, (0.0, 1.0))
    @test_throws ["single unknown"] IntervalNonlinearFunction(sys)
    @test_throws ["single unknown"] IntervalNonlinearProblem(
        sys, (0.0, 1.0); expression = Val{true})
    @test_throws ["single unknown"] IntervalNonlinearFunction(
        sys; expression = Val{true})
end

@testset "Vector parameter used unscalarized and partially scalarized" begin
    @variables x y
    @parameters p[1:2] (f::Function)(..)

    @mtkcompile sys = System([x^2 - p[1]^2 ~ 0, y^2 ~ f(p)])
    @test !any(isequal(p[1]), parameters(sys))
    @test is_parameter(sys, p)
end

@testset "Can convert from `System`" begin
    @variables x(t) y(t)
    @parameters p q r
    @named sys = System([D(x) ~ p * x^3 + q, 0 ~ -y + q * x - r, r ~ 3p], t;
        defaults = [x => 1.0, p => missing], guesses = [p => 1.0],
        initialization_eqs = [p^3 + q^3 ~ 4r])
    nlsys = NonlinearSystem(sys)
    nlsys = complete(nlsys)
    defs = defaults(nlsys)
    @test length(defs) == 6
    @test defs[x] == 1.0
    @test defs[p] === missing
    @test isinf(defs[t])
    @test length(guesses(nlsys)) == 1
    @test guesses(nlsys)[p] == 1.0
    @test length(initialization_equations(nlsys)) == 1
    @test length(parameter_dependencies(nlsys)) == 1
    @test length(equations(nlsys)) == 2
    @test all(iszero, [eq.lhs for eq in equations(nlsys)])
    @test nameof(nlsys) == nameof(sys)
    @test ModelingToolkit.iscomplete(nlsys)

    sys1 = complete(sys; split = false)
    nlsys = NonlinearSystem(sys1)
    @test ModelingToolkit.iscomplete(nlsys)
    @test !ModelingToolkit.is_split(nlsys)

    sys2 = complete(sys)
    nlsys = NonlinearSystem(sys2)
    @test ModelingToolkit.iscomplete(nlsys)
    @test ModelingToolkit.is_split(nlsys)

    sys3 = mtkcompile(sys)
    nlsys = NonlinearSystem(sys3)
    @test length(equations(nlsys)) == length(ModelingToolkit.observed(nlsys)) == 1

    prob = NonlinearProblem(sys3, [q => 2.0])
    @test prob.f.initialization_data.initializeprobmap === nothing
    sol = solve(prob)
    @test SciMLBase.successful_retcode(sol)
    @test sol.ps[p ^ 3 + q ^ 3]≈sol.ps[4r] atol=1e-10

    @testset "Differential inside expression also substituted" begin
        @named sys = System([0 ~ y * D(x) + x^2 - p, 0 ~ x * D(y) + y * p], t)
        nlsys = NonlinearSystem(sys)
        vs = ModelingToolkit.vars(equations(nlsys))
        @test !in(D(x), vs)
        @test !in(D(y), vs)
    end
end

@testset "oop `NonlinearLeastSquaresProblem` with `u0 === nothing`" begin
    @variables x y
    @named sys = System([0 ~ x - y], [], []; observed = [x ~ 1.0, y ~ 1.0])
    prob = NonlinearLeastSquaresProblem{false}(complete(sys), nothing)
    sol = solve(prob)
    resid = sol.resid
    @test resid == [0.0]
    @test resid isa Vector
    prob = NonlinearLeastSquaresProblem{false}(
        complete(sys), nothing; u0_constructor = splat(SVector))
    sol = solve(prob)
    resid = sol.resid
    @test resid == [0.0]
    @test resid isa SVector
end

@testset "`ProblemTypeCtx`" begin
    @variables x
    @mtkcompile sys = System(
        [0 ~ x^2 - 4x + 4]; metadata = [ModelingToolkit.ProblemTypeCtx => "A"])
    prob = NonlinearProblem(sys, [x => 1.0])
    @test prob.problem_type == "A"
end
