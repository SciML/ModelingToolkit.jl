using ModelingToolkit, StaticArrays, LinearAlgebra
using ModelingToolkit: get_metadata
using DiffEqBase, SparseArrays
using Test
using NonlinearSolve
using ModelingToolkit: value
using ModelingToolkit: get_default_or_guess

canonequal(a, b) = isequal(simplify(a), simplify(b))

# Define some variables
@parameters t σ ρ β
@constants h = 1
@variables x y z

function test_nlsys_inference(name, sys, vs, ps)
    @testset "NonlinearSystem construction: $name" begin
        @test Set(states(sys)) == Set(value.(vs))
        @test Set(parameters(sys)) == Set(value.(ps))
    end
end

# Define a nonlinear system
eqs = [0 ~ σ * (y - x) * h,
    0 ~ x * (ρ - z) - y,
    0 ~ x * y - β * z]
@named ns = NonlinearSystem(eqs, [x, y, z], [σ, ρ, β], defaults = Dict(x => 2))
@test eval(toexpr(ns)) == ns
test_nlsys_inference("standard", ns, (x, y, z), (σ, ρ, β))
@test begin
    f = eval(generate_function(ns, [x, y, z], [σ, ρ, β])[2])
    du = [0.0, 0.0, 0.0]
    f(du, [1, 2, 3], [1, 2, 3])
    du ≈ [1, -3, -7]
end

# Now nonlinear system with only variables
@variables x y z
@parameters σ ρ β

# Define a nonlinear system
eqs = [0 ~ σ * (y - x),
    y ~ x * (ρ - z),
    β * z ~ x * y]
@named ns = NonlinearSystem(eqs, [x, y, z], [σ, ρ, β])
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
nlsys_func = generate_function(ns, [x, y, z], [σ, ρ, β])
jac_func = generate_jacobian(ns)
f = @eval eval(nlsys_func)

# Intermediate calculations
a = y - x
# Define a nonlinear system
eqs = [0 ~ σ * a * h,
    0 ~ x * (ρ - z) - y,
    0 ~ x * y - β * z]
@named ns = NonlinearSystem(eqs, [x, y, z], [σ, ρ, β])
nlsys_func = generate_function(ns, [x, y, z], [σ, ρ, β])
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

prob = NonlinearProblem(ns, ones(3), ones(3))
@test prob.f.sys === ns
sol = solve(prob, NewtonRaphson())
@test sol.u[1] ≈ sol.u[2]

@test_throws ArgumentError NonlinearProblem(ns, ones(4), ones(3))

@variables u F s a
eqs1 = [
    0 ~ σ * (y - x) * h + F,
    0 ~ x * (ρ - z) - u,
    0 ~ x * y - β * z,
    0 ~ x + y - z - u,
]

lorenz = name -> NonlinearSystem(eqs1, [x, y, z, u, F], [σ, ρ, β], name = name)
lorenz1 = lorenz(:lorenz1)
@test_throws ArgumentError NonlinearProblem(lorenz1, zeros(5))
lorenz2 = lorenz(:lorenz2)
@named connected = NonlinearSystem([s ~ a + lorenz1.x
        lorenz2.y ~ s * h
        lorenz1.F ~ lorenz2.u
        lorenz2.F ~ lorenz1.u], [s, a], [],
    systems = [lorenz1, lorenz2])
@test_nowarn alias_elimination(connected)

# system promotion
using OrdinaryDiffEq
@variables t
D = Differential(t)
@named subsys = convert_system(ODESystem, lorenz1, t)
@named sys = ODESystem([D(subsys.x) ~ subsys.x + subsys.x], t, systems = [subsys])
sys = structural_simplify(sys)
u0 = [subsys.x => 1, subsys.z => 2.0, subsys.y => 1.0]
prob = ODEProblem(sys, u0, (0, 1.0), [subsys.σ => 1, subsys.ρ => 2, subsys.β => 3])
sol = solve(prob, FBDF(), reltol = 1e-7, abstol = 1e-7)
@test sol[subsys.x] + sol[subsys.y] - sol[subsys.z]≈sol[subsys.u] atol=1e-7
@test_throws ArgumentError convert_system(ODESystem, sys, t)

@parameters t σ ρ β
@variables x y z

# Define a nonlinear system
eqs = [0 ~ σ * (y - x),
    0 ~ x * (ρ - z) - y,
    0 ~ x * y - β * z * h]
@named ns = NonlinearSystem(eqs, [x, y, z], [σ, ρ, β])
np = NonlinearProblem(ns, [0, 0, 0], [1, 2, 3], jac = true, sparse = true)
@test calculate_jacobian(ns, sparse = true) isa SparseMatrixCSC

# issue #819
@testset "Combined system name collisions" begin
    function makesys(name)
        @parameters a
        @variables x f

        NonlinearSystem([0 ~ -a * x + f], [x, f], [a]; name)
    end

    function issue819()
        sys1 = makesys(:sys1)
        sys2 = makesys(:sys1)
        @test_throws ArgumentError NonlinearSystem([sys2.f ~ sys1.x, sys1.f ~ 0], [], [],
            systems = [sys1, sys2], name = :foo)
    end
    issue819()
end

# issue #1115
@testset "Extending a NonlinearSystem with no iv" begin
    @parameters a b
    @variables x y
    eqs1 = [
        0 ~ a * x,
    ]
    eqs2 = [
        0 ~ b * y,
    ]

    @named sys1 = NonlinearSystem(eqs1, [x], [a])
    @named sys2 = NonlinearSystem(eqs2, [y], [b])
    @named sys3 = extend(sys1, sys2)

    @test isequal(union(Set(parameters(sys1)), Set(parameters(sys2))),
        Set(parameters(sys3)))
    @test isequal(union(Set(states(sys1)), Set(states(sys2))), Set(states(sys3)))
    @test isequal(union(Set(equations(sys1)), Set(equations(sys2))), Set(equations(sys3)))
end

# observed variable handling
@variables t x(t) RHS(t)
@parameters τ
@named fol = NonlinearSystem([0 ~ (1 - x * h) / τ], [x], [τ];
    observed = [RHS ~ (1 - x) / τ])
@test isequal(RHS, @nonamespace fol.RHS)
RHS2 = RHS
@unpack RHS = fol
@test isequal(RHS, RHS2)

# issue #1358
@variables t
@variables v1(t) v2(t) i1(t) i2(t)
eq = [v1 ~ sin(2pi * t * h)
    v1 - v2 ~ i1
    v2 ~ i2
    i1 ~ i2]
@named sys = ODESystem(eq)
@test length(equations(structural_simplify(sys))) == 0

@testset "Issue: 1504" begin
    @variables u[1:4]

    eqs = [u[1] ~ 1,
        u[2] ~ 1,
        u[3] ~ 1,
        u[4] ~ h]

    sys = NonlinearSystem(eqs, collect(u[1:4]), Num[], defaults = Dict([]), name = :test)
    prob = NonlinearProblem(sys, ones(length(states(sys))))

    sol = NonlinearSolve.solve(prob, NewtonRaphson())

    @test sol[u] ≈ ones(4)
end

@variables x(t)
@parameters a
eqs = [0 ~ a * x]

testdict = Dict([:test => 1])
@named sys = NonlinearSystem(eqs, [x], [a], metadata = testdict)
@test get_metadata(sys) == testdict

@testset "Remake" begin
    @parameters a=1.0 b=1.0 c=1.0
    @constants h = 1
    @variables x y z

    eqs = [0 ~ a * (y - x) * h,
        0 ~ x * (b - z) - y,
        0 ~ x * y - c * z]
    @named sys = NonlinearSystem(eqs, [x, y, z], [a, b, c], defaults = Dict(x => 2.0))
    prob = NonlinearProblem(sys, ones(length(states(sys))))

    prob_ = remake(prob, u0 = [1.0, 2.0, 3.0], p = [1.1, 1.2, 1.3])
    @test prob_.u0 == [1.0, 2.0, 3.0]
    @test prob_.p == [1.1, 1.2, 1.3]

    prob_ = remake(prob, u0 = Dict(y => 2.0), p = Dict(a => 2.0))
    @test prob_.u0 == [1.0, 2.0, 1.0]
    @test prob_.p == [2.0, 1.0, 1.0]
end

@testset "Initialization System" begin
    # Define the Lotka Volterra system which begins at steady state
    @parameters t
    pars = @parameters a=1.5 b=1.0 c=3.0 d=1.0 dx_ss=1e-5

    vars = @variables begin
        dx(t),
        dy(t),
        (x(t) = dx ~ dx_ss), [guess = 0.5]
        (y(t) = dy ~ 0), [guess = -0.5]
    end

    D = Differential(t)

    eqs = [dx ~ a * x - b * x * y
        dy ~ -c * y + d * x * y
        D(x) ~ dx
        D(y) ~ dy]

    @named sys = ODESystem(eqs, t, vars, pars)

    sys_simple = structural_simplify(sys)

    # Set up the initialization system
    sys_init = initializesystem(sys_simple)

    sys_init_simple = structural_simplify(sys_init)

    prob = NonlinearProblem(sys_init_simple, get_default_or_guess.(states(sys_init_simple)))

    @test prob.u0 == [0.5, -0.5]

    sol = solve(prob)
    @test sol.retcode == SciMLBase.ReturnCode.Success

    # Confirm for all the states of the non-simplified system
    @test all(.≈(sol[states(sys)], [1e-5, 0, 1e-5 / 1.5, 0]; atol = 1e-8))

    # Confirm for all the states of the simplified system
    @test all(.≈(sol[states(sys_simple)], [1e-5 / 1.5, 0]; atol = 1e-8))
end
